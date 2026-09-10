exposure_labels <- function(windowed = FALSE) {
  if (!windowed) return(c("other_observed_time", "pregnancy", "postpartum_days_1_180"))
  c("other_observed_time", "pregnancy", paste0("postpartum_days_", seq(1, 151, 30), "_", seq(30, 180, 30)))
}

model_matrix <- function(rows, variables, windowed = FALSE, intercept = TRUE, drop_absent = intercept, labels = exposure_labels(windowed)) {
  definitions <- data.table::data.table(name = labels[-1], column = "risk_code", code = seq_along(labels)[-1] - 1L)
  add <- function(nms, col, codes) data.table::data.table(name = nms, column = col, code = codes)
  specs <- list(
    age = add(c("age_lt20", "age_30_39", "age_40_49"), "age_code", 1:3),
    calendar = add(c("calendar_2015_2019", "calendar_2020_2025"), "calendar_code", 1:2),
    calendar_utt = add("calendar_2020_2025", "calendar_code", 1L),
    who = add(c("who_stage_3_4", "who_stage_missing"), "who_code", 1:2),
    art_duration = add(c("art_duration_lt6m", "art_duration_6_11m", "art_duration_1_1_99y", "art_duration_2_4_99y"), "art_duration_code", 0:3),
    cd4 = add(c("cd4_350_499", "cd4_200_349", "cd4_lt200", "cd4_missing"), "cd4_code", 1:4),
    bmi = add(c("bmi_underweight", "bmi_ge25", "bmi_missing"), "bmi_code", 1:3),
    marital = add(c("married_or_partnered", "marital_missing"), "marital_code", 1:2),
    tpt = add("recorded_tpt_before_entry", "tpt_before_entry", 1L))
  if ("clinic" %in% variables) {
    counts <- rows[, .N, by = clinic_code]
    reference <- counts$clinic_code[which.max(counts$N)]
    codes <- setdiff(sort(counts$clinic_code), reference)
    specs$clinic <- add(paste0("clinic_", codes), "clinic_code", codes)
  }
  if (any(!variables %in% names(specs))) stop("Unknown model covariate.")
  definitions <- data.table::rbindlist(c(list(definitions), specs[variables]))
  x <- vapply(seq_len(nrow(definitions)), function(j)
    as.numeric(rows[[definitions$column[j]]] == definitions$code[j]), numeric(nrow(rows)))
  colnames(x) <- definitions$name
  if (anyNA(x)) stop("Model matrix contains missing covariates.")
  if (drop_absent) x <- x[, colSums(abs(x)) > 0, drop = FALSE]
  if (intercept) x <- cbind(intercept = 1, x)
  x
}

package_model <- function(beta, covariance, information, score, rows, method, converged, iterations, loglik) {
  se <- sqrt(pmax(diag(covariance), 0))
  estimates <- data.table::data.table(feature = names(beta), log_irr = as.numeric(beta), standard_error = se,
    irr = exp(beta), ci_lower = exp(beta - 1.96 * se), ci_upper = exp(beta + 1.96 * se),
    p_value = 2 * stats::pnorm(-abs(beta / se)))
  rank <- qr(information, tol = 1e-10)$rank
  if (!converged || rank != length(beta) || any(!is.finite(covariance)) || max(abs(score)) > 1e-5)
    stop("Model failed convergence, rank, score, or covariance checks: ", method)
  list(beta = beta, covariance = covariance, estimates = estimates,
    diagnostics = list(method = method, converged = converged, iterations = iterations,
      information_rank = rank, parameter_count = length(beta), maximum_absolute_score = max(abs(score)),
      analysis_n = data.table::uniqueN(rows$patient_id), rows = nrow(rows), events = sum(rows$events),
      person_days = sum(rows$days), log_likelihood = loglik))
}

fit_cohort_poisson <- function(rows, variables = c("age", "calendar"), windowed = FALSE,
                               labels = exposure_labels(windowed), x = NULL) {
  stopifnot(all(rows$days > 0), all(rows$events >= 0))
  if (is.null(x)) x <- model_matrix(rows, variables, windowed, labels = labels)
  y <- rows$events; days <- rows$days
  start <- rep(0, ncol(x)); start[1] <- log(sum(y) / sum(days))
  fit <- stats::glm.fit(x, y, family = stats::poisson(), offset = log(days), start = start,
    control = stats::glm.control(epsilon = 1e-12, maxit = 100))
  beta <- fit$coefficients; mu <- fit$fitted.values
  information <- crossprod(x, x * mu); bread <- solve(information)
  cluster_scores <- rowsum(x * (y - mu), rows$patient_id, reorder = FALSE)
  G <- nrow(cluster_scores); N <- nrow(x); P <- ncol(x)
  stopifnot(G > 1L, N > P)
  correction <- G / (G - 1) * (N - 1) / (N - P)
  covariance <- correction * bread %*% crossprod(cluster_scores) %*% bread
  result <- package_model(beta, covariance, information, crossprod(x, y - mu), rows,
    "R stats::glm.fit Poisson; participant-clustered sandwich with original finite-sample correction",
    fit$converged, fit$iter, sum(y * (log(days) + x %*% beta) - mu))
  result$diagnostics$finite_sample_correction <- correction
  result
}

conditional_poisson_state <- function(beta, x, days, events, groups) {
  eta <- as.numeric(x %*% beta) + log(days)
  maxima <- as.numeric(tapply(eta, groups, max))
  weights <- exp(eta - maxima[groups])
  totals <- as.numeric(rowsum(weights, groups, reorder = TRUE))
  probabilities <- weights / totals[groups]
  means <- rowsum(x * probabilities, groups, reorder = TRUE)
  score <- colSums(x * events) - colSums(means)
  information <- crossprod(x, x * probabilities) - crossprod(means)
  # Omit the event log(days) term; it is constant with respect to beta.
  loglik <- sum(events * as.numeric(x %*% beta)) - sum(maxima + log(totals))
  list(loglik = loglik, score = score, information = information)
}

fit_sccs_poisson <- function(rows, variables = c("age", "calendar"), windowed = FALSE,
                             labels = exposure_labels(windowed), x = NULL) {
  stopifnot(all(rows$days > 0), all(rows[, sum(events), by = patient_id]$V1 == 1))
  if (is.null(x)) x <- model_matrix(rows, variables, windowed, intercept = FALSE, drop_absent = FALSE, labels = labels)
  groups <- match(rows$patient_id, unique(rows$patient_id))
  state <- function(beta) conditional_poisson_state(beta, x, rows$days, rows$events, groups)
  beta <- setNames(rep(0, ncol(x)), colnames(x)); converged <- FALSE
  for (iteration in seq_len(100L)) {
    s <- state(beta); step <- solve(s$information, s$score)
    if (max(abs(step)) < 1e-10) { converged <- TRUE; break }
    scale <- 1
    while (state(beta + scale * step)$loglik < s$loglik - 1e-10 && scale > 1e-8) scale <- scale / 2
    beta <- beta + scale * step
    if (max(abs(scale * step)) < 1e-9) { converged <- TRUE; break }
  }
  s <- state(beta)
  package_model(beta, solve(s$information), s$information, s$score, rows,
    "Participant-conditioned Poisson; inverse-information covariance", converged, iteration, s$loglik)
}
