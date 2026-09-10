# Postpartum restricted cubic spline, centered at day 90.
postpartum_spline_basis <- function(days, centered = TRUE) {
  basis <- function(x) {
    knots <- c(10, 64, 117, 171); last <- knots[4]; penultimate <- knots[3]
    nonlinear <- matrix(vapply(knots[1:2], function(knot) {
      (pmax(x - knot, 0)^3 - pmax(x - penultimate, 0)^3 * (last - knot) / (last - penultimate) +
        pmax(x - last, 0)^3 * (penultimate - knot) / (last - penultimate)) / ((last - knots[1])^2 * 30)
    }, numeric(length(x))), nrow = length(x), ncol = 2L)
    cbind(postpartum_rcs_linear = x / 30, postpartum_rcs_nonlinear_1 = nonlinear[, 1], postpartum_rcs_nonlinear_2 = nonlinear[, 2])
  }
  result <- basis(days)
  if (centered) result <- sweep(result, 2, basis(90)[1, ], "-")
  result
}

spline_person_time <- function(intervals, episodes, sccs = FALSE) {
  x <- data.table::copy(data.table::as.data.table(intervals)[included == TRUE])
  if (sccs) x <- x[primary_case == TRUE]
  other <- x[risk_code != 2L, .(patient_id, state = risk_code, age_code, calendar_code, days, events)]
  pp <- x[risk_code == 2L]
  e <- data.table::as.data.table(episodes)
  idx <- match(paste(pp$patient_id, pp$documented_pregnancy_order, sep = ":"),
    paste(e$patient_id, e$documented_pregnancy_order, sep = ":"))
  first_day <- as.integer(pp$interval_start - e$end[idx])
  stopifnot(all(first_day >= 1L), all(first_day + pp$days - 1L <= 180L), all(pp$events == 0L | pp$days == 1L))
  rows <- rep(seq_len(nrow(pp)), pp$days)
  detailed <- data.table::data.table(patient_id = pp$patient_id[rows], state = first_day[rows] + sequence(pp$days),
    age_code = pp$age_code[rows], calendar_code = pp$calendar_code[rows], days = 1L, events = pp$events[rows])
  result <- data.table::rbindlist(list(other, detailed))[, .(days = sum(days), events = sum(events)),
    by = .(patient_id, state, age_code, calendar_code)]
  data.table::setorder(result, patient_id, state, age_code, calendar_code)
  result
}

spline_design_matrix <- function(rows, intercept = TRUE) {
  postpartum <- rows$state >= 2L
  b <- postpartum_spline_basis(pmax(rows$state - 1L, 1L))
  b[!postpartum, ] <- 0
  x <- cbind(pregnancy = as.numeric(rows$state == 1L), postpartum_at_day_90 = as.numeric(postpartum), b,
    age_lt20 = as.numeric(rows$age_code == 1L), age_30_39 = as.numeric(rows$age_code == 2L),
    age_40_49 = as.numeric(rows$age_code == 3L), calendar_2015_2019 = as.numeric(rows$calendar_code == 1L),
    calendar_2020_2025 = as.numeric(rows$calendar_code == 2L))
  if (intercept) x <- cbind(intercept = 1, x)
  x
}

postpartum_spline_predictions <- function(model) {
  days <- seq_len(180L)
  x <- matrix(0, nrow = length(days), ncol = length(model$beta), dimnames = list(NULL, names(model$beta)))
  x[, "postpartum_at_day_90"] <- 1
  x[, c("postpartum_rcs_linear", "postpartum_rcs_nonlinear_1", "postpartum_rcs_nonlinear_2")] <- postpartum_spline_basis(days)
  log_irr <- as.numeric(x %*% model$beta)
  se <- sqrt(rowSums((x %*% model$covariance) * x))
  data.table::data.table(postpartum_day = days, log_irr = log_irr, standard_error = se,
    irr = exp(log_irr), ci_lower = exp(log_irr - 1.96 * se), ci_upper = exp(log_irr + 1.96 * se))
}

postpartum_spline_tests <- function(model) {
  names <- c("postpartum_at_day_90", "postpartum_rcs_linear", "postpartum_rcs_nonlinear_1", "postpartum_rcs_nonlinear_2")
  groups <- list(overall_postpartum_curve = names, variation_over_postpartum_time = names[-1], nonlinearity = names[3:4])
  data.table::rbindlist(lapply(names(groups), function(label) {
    selection <- diag(length(model$beta))[match(groups[[label]], names(model$beta)), , drop = FALSE]
    wald_contrast(model, selection, label)
  }))
}
