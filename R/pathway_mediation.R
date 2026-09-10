# Contact mediation with woman-clustered joint covariance.
mediation_design <- function(rows, specification = "age_calendar") {
  x <- cbind(intercept = 1, pregnancy = as.numeric(rows$risk_code == 1L),
    postpartum_days_1_180 = as.numeric(rows$risk_code == 2L), age_lt20 = as.numeric(rows$age_code == 1L),
    age_30_39 = as.numeric(rows$age_code == 2L), age_40_49 = as.numeric(rows$age_code == 3L),
    calendar_2015_2019 = as.numeric(rows$calendar_code == 1L), calendar_2020_2025 = as.numeric(rows$calendar_code == 2L))
  if (specification == "expanded") {
    art <- vapply(0:3, function(code) as.numeric(rows$art_duration_code == code), numeric(nrow(rows)))
    colnames(art) <- c("art_lt6m", "art_6_11m", "art_1_1_99y", "art_2_4_99y")
    clinic_days <- rows[, .(days = sum(days)), by = clinic_code][order(clinic_code)]
    reference <- clinic_days$clinic_code[which.max(clinic_days$days)]
    clinic_codes <- setdiff(sort(unique(rows$clinic_code)), reference)
    clinic <- vapply(clinic_codes, function(code) as.numeric(rows$clinic_code == code), numeric(nrow(rows)))
    colnames(clinic) <- paste0("clinic_", clinic_codes)
    x <- cbind(x, art, clinic, partnered = as.numeric(rows$marital_code == 1L),
      partner_missing = as.numeric(rows$marital_code == 2L), tpt_before_entry = rows$tpt_before_entry,
      cd4_350_499 = as.numeric(rows$cd4_code == 1L), cd4_200_349 = as.numeric(rows$cd4_code == 2L),
      cd4_lt200 = as.numeric(rows$cd4_code == 3L), cd4_missing = as.numeric(rows$cd4_code == 4L),
      bmi_underweight = as.numeric(rows$bmi_code == 1L), bmi_ge25 = as.numeric(rows$bmi_code == 2L),
      bmi_missing = as.numeric(rows$bmi_code == 3L), pre_entry_visit_rate_log = rows$pre_entry_visit_rate_log)
  } else if (specification != "age_calendar") stop("Unknown mediation covariate specification.")
  stopifnot(!anyNA(x))
  x
}

mediation_component <- function(x, y, days, patient_id, family) {
  if (family == "linear") {
    information <- crossprod(x, x * days)
    beta <- as.numeric(solve(information, crossprod(x, days * y)))
    names(beta) <- colnames(x)
    fitted <- as.numeric(x %*% beta)
    residual <- days * (y - fitted)
    converged <- TRUE; iterations <- 1L
    loglik <- NA_real_
  } else {
    start <- rep(0, ncol(x)); start[1] <- log(max(sum(y), 0.5) / sum(days))
    fit <- stats::glm.fit(x, y, family = stats::poisson(), offset = log(days), start = start,
      control = stats::glm.control(epsilon = 1e-12, maxit = 100L))
    beta <- fit$coefficients; fitted <- fit$fitted.values
    information <- crossprod(x, x * fitted)
    residual <- y - fitted; converged <- fit$converged; iterations <- fit$iter
    loglik <- sum(y * (log(days) + as.numeric(x %*% beta)) - fitted)
  }
  bread <- solve(information)
  scores <- rowsum(x * residual, patient_id, reorder = FALSE)
  influence <- scores %*% t(bread)
  G <- nrow(scores); correction <- G / (G - 1)
  covariance <- correction * crossprod(influence)
  dimnames(covariance) <- list(colnames(x), colnames(x))
  score <- colSums(x * residual)
  rank <- qr(information, tol = 1e-10)$rank
  if (!converged || rank != ncol(x) || any(!is.finite(covariance)) || max(abs(score)) > 1e-5)
    stop("Mediation component failed convergence, rank, score, or covariance checks: ", family)
  list(beta = beta, covariance = covariance, influence = influence,
    diagnostics = list(family = family, converged = converged, iterations = iterations,
      information_rank = rank, parameters = ncol(x), rows = nrow(x), patients = G,
      maximum_absolute_score = max(abs(score)), finite_sample_correction = correction, log_likelihood = loglik))
}

mediation_scenario <- function(x, days, alpha, beta, outcome_exposure, mediator_exposure) {
  xm <- xo <- x
  xm[, "pregnancy"] <- as.numeric(mediator_exposure == 1L)
  xm[, "postpartum_days_1_180"] <- as.numeric(mediator_exposure == 2L)
  linear <- as.numeric(xm %*% alpha)
  active <- linear > 0
  mediator <- pmax(linear, 0)
  xo[, "pregnancy"] <- as.numeric(outcome_exposure == 1L)
  xo[, "postpartum_days_1_180"] <- as.numeric(outcome_exposure == 2L)
  z <- cbind(xo, mediator, pregnancy_mediator = xo[, "pregnancy"] * mediator,
    postpartum_mediator = xo[, "postpartum_days_1_180"] * mediator)
  rate <- exp(pmin(pmax(as.numeric(z %*% beta), -35), 25))
  weighted <- days * rate; normalized <- weighted / sum(weighted)
  slope <- beta[ncol(x) + 1L] + if (outcome_exposure == 0L) 0 else beta[ncol(x) + 1L + outcome_exposure]
  # Predictions below zero are truncated on the log(1 + rate) scale.
  gradient_alpha <- as.numeric(slope * crossprod(normalized, xm * active))
  gradient_beta <- as.numeric(crossprod(normalized, z))
  list(log_rate = log(sum(weighted) / sum(days)), gradient = c(gradient_alpha, gradient_beta),
    truncated_mediator_rows = sum(!active))
}

mediation_delta <- function(log_effect, gradient, covariance) {
  se <- sqrt(max(as.numeric(crossprod(gradient, covariance %*% gradient)), 0))
  list(rr = exp(log_effect), ci_lower = exp(log_effect - 1.96 * se), ci_upper = exp(log_effect + 1.96 * se),
    log_rr = log_effect, log_rr_standard_error = se)
}

fit_contact_mediation <- function(rows, mediator_count_column, specification = "age_calendar") {
  x <- mediation_design(rows, specification)
  days <- rows$days; event <- rows$events
  mediator <- log1p(rows[[mediator_count_column]] / days * 100)
  alpha <- mediation_component(x, mediator, days, rows$patient_id, "linear")
  total <- mediation_component(x, event, days, rows$patient_id, "poisson")
  z <- cbind(x, mediator, pregnancy_mediator = x[, "pregnancy"] * mediator,
    postpartum_mediator = x[, "postpartum_days_1_180"] * mediator)
  outcome <- mediation_component(z, event, days, rows$patient_id, "poisson")
  interaction <- tail(seq_len(ncol(z)), 2L)
  ib <- outcome$beta[interaction]; iv <- outcome$covariance[interaction, interaction]
  interaction_wald <- as.numeric(crossprod(ib, solve(iv, ib)))
  G <- data.table::uniqueN(rows$patient_id)
  cross <- G / (G - 1) * crossprod(alpha$influence, outcome$influence)
  joint <- rbind(cbind(alpha$covariance, cross), cbind(t(cross), outcome$covariance))
  joint_names <- c(paste0("mediator:", colnames(x)), paste0("outcome:", colnames(z)))
  dimnames(joint) <- list(joint_names, joint_names)
  baseline <- mediation_scenario(x, days, alpha$beta, outcome$beta, 0L, 0L)
  scenarios <- list(reference = baseline)
  effects <- lapply(1:2, function(code) {
    label <- c("pregnancy", "postpartum_days_1_180")[code]
    a0 <- mediation_scenario(x, days, alpha$beta, outcome$beta, code, 0L)
    aa <- mediation_scenario(x, days, alpha$beta, outcome$beta, code, code)
    scenarios[[paste0(label, "_reference_mediator")]] <<- a0
    scenarios[[paste0(label, "_exposure_mediator")]] <<- aa
    direct <- mediation_delta(a0$log_rate - baseline$log_rate, a0$gradient - baseline$gradient, joint)
    indirect <- mediation_delta(aa$log_rate - a0$log_rate, aa$gradient - a0$gradient, joint)
    decomposed <- mediation_delta(aa$log_rate - baseline$log_rate, aa$gradient - baseline$gradient, joint)
    tlog <- total$beta[label]; tse <- sqrt(total$covariance[label, label])
    result <- list(exposure = label, mediator_count_definition = mediator_count_column, covariate_specification = specification,
      total_effect_rr = exp(tlog), total_effect_ci_lower = exp(tlog - 1.96 * tse), total_effect_ci_upper = exp(tlog + 1.96 * tse),
      direct_effect_rr = direct$rr, direct_effect_ci_lower = direct$ci_lower, direct_effect_ci_upper = direct$ci_upper,
      indirect_effect_rr = indirect$rr, indirect_effect_ci_lower = indirect$ci_lower, indirect_effect_ci_upper = indirect$ci_upper,
      decomposed_total_rr = decomposed$rr, decomposed_total_ci_lower = decomposed$ci_lower, decomposed_total_ci_upper = decomposed$ci_upper,
      exposure_mediator_coefficient = alpha$beta[label], baseline_mediator_outcome_log_rr_coefficient = outcome$beta[ncol(x) + 1L],
      exposure_mediator_interaction_coefficient = outcome$beta[ncol(x) + 1L + code],
      joint_exposure_mediator_interaction_p = stats::pchisq(interaction_wald, 2, lower.tail = FALSE),
      patients = G, events = sum(event), person_days = sum(days), mediator_model_converged = TRUE,
      total_model_converged = total$diagnostics$converged, outcome_model_converged = outcome$diagnostics$converged,
      interaction_model_converged = outcome$diagnostics$converged)
    data.table::as.data.table(result)
  })
  list(effects = data.table::rbindlist(effects), mediator = alpha, total = total, outcome = outcome,
    joint_covariance = joint, scenarios = scenarios, interaction_wald = interaction_wald,
    specification = list(mediator = mediator_count_column, covariates = specification,
      standardization = "Observed person-day distribution", interval_method = "Joint cluster-robust multivariate delta method"))
}
