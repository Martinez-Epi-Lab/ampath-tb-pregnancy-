# Secondary exposure definitions and sensitivity analyses.
secondary_labels <- function(specification) {
  switch(specification,
    order = c("other_observed_time", "pregnancy_first_documented", "postpartum_1_180_first_documented",
      "pregnancy_second_documented", "postpartum_1_180_second_documented",
      "pregnancy_third_or_later_documented", "postpartum_1_180_third_or_later_documented"),
    art_timing = c("other_observed_time", "pregnancy_underway_at_art", "postpartum_1_180_underway_at_art",
      "pregnancy_conceived_after_art", "postpartum_1_180_conceived_after_art"),
    fully_observed = c("unexposed", "pregnancy", "postpartum", "pregnancy_nonqualifying_nuisance", "postpartum_nonqualifying_nuisance"),
    live_birth = c("unexposed", "pregnancy_live_birth", "postpartum_live_birth", "pregnancy_other_outcome_nuisance", "postpartum_other_outcome_nuisance"),
    postpartum_31_180 = c("other_observed_time", "pregnancy", "postpartum_days_1_30", "postpartum_days_31_180"),
    postpartum_60day = c("other_observed_time", "pregnancy", "postpartum_days_1_60", "postpartum_days_61_120", "postpartum_days_121_180"),
    exposure_labels())
}

split_intervals_at_dates <- function(intervals, cuts) {
  x <- data.table::copy(data.table::as.data.table(intervals))
  if (!nrow(cuts) || !nrow(x)) return(x)
  x[, interval_id := seq_len(.N)]
  cut_table <- unique(data.table::as.data.table(cuts)[!is.na(cut_date), .(patient_id, cut_date)])
  extra <- x[cut_table, on = .(patient_id, interval_start < cut_date, interval_end >= cut_date),
    nomatch = 0, allow.cartesian = TRUE, .(interval_id, left = i.cut_date)]
  starts <- unique(data.table::rbindlist(list(x[, .(interval_id, left = interval_start)], extra)))
  data.table::setorder(starts, interval_id, left)
  starts[, right := data.table::shift(left, type = "lead") - 1L, by = interval_id]
  starts[is.na(right), right := x$interval_end[interval_id]]
  out <- x[starts$interval_id]
  out[, events := as.integer(events == 1L & starts$left == interval_start)]
  out[, `:=`(interval_start = starts$left, interval_end = starts$right,
    days = as.integer(starts$right - starts$left) + 1L, interval_id = NULL)]
  out
}

secondary_rows <- function(intervals, cohort, episodes, specification = "primary", sccs = FALSE, outcomes = NULL, art_duration = FALSE) {
  x <- data.table::copy(data.table::as.data.table(intervals))
  p <- data.table::as.data.table(episodes)
  all_orders <- specification %in% c("all_pregnancies", "order", "art_timing")
  x[, included := all_orders | !(order_code %in% c(5L, 6L))]
  if (specification == "exclude_preconception") {
    pre <- p[documented_pregnancy_order <= 2L, .(patient_id, pre_start = original_start - 180L, pre_end = original_start - 1L)]
    cuts <- data.table::rbindlist(list(pre[, .(patient_id, cut_date = pre_start)], pre[, .(patient_id, cut_date = pre_end + 1L)]))
    x <- split_intervals_at_dates(x, cuts)
    x[, row_id := seq_len(.N)]
    remove <- x[pre, on = .(patient_id, interval_start >= pre_start, interval_end <= pre_end), nomatch = 0, allow.cartesian = TRUE, row_id]
    x[row_id %in% remove & risk_code == 0L, included := FALSE]
    x[, row_id := NULL]
  }
  if (specification %in% c("utt_era", "cutoff_2024")) {
    boundary <- as.Date(if (specification == "utt_era") "2016-07-01" else "2024-12-31")
    if (specification == "utt_era") {
      x <- x[interval_end >= boundary]
      x[, interval_start := pmax(interval_start, boundary)]
      x[, calendar_code := as.integer(interval_start >= as.Date("2020-01-01"))]
    } else {
      x <- x[interval_start <= boundary]
      x[, interval_end := pmin(interval_end, boundary)]
    }
    x[, days := as.integer(interval_end - interval_start) + 1L]
  }
  if (specification == "order") x[, risk_code := order_code]
  if (specification %in% c("art_timing", "fully_observed", "live_birth")) {
    keys <- paste(p$patient_id, p$documented_pregnancy_order, sep = ":")
    idx <- match(paste(x$patient_id, x$documented_pregnancy_order, sep = ":"), keys)
    if (specification == "art_timing") {
      after <- p$start > cohort$arvstartTHERAP[match(p$patient_id, cohort$patient_id)]
      x[which(x$risk_code > 0 & after[idx]), risk_code := risk_code + 2L]
    }
    if (specification == "fully_observed") {
      ci <- match(p$patient_id, cohort$patient_id)
      qualifies <- p$start >= cohort$entry[ci] & p$end <= cohort$administrative_exit[ci]
      x[which(x$risk_code > 0 & !qualifies[idx]), risk_code := risk_code + 2L]
    }
    if (specification == "live_birth") {
      stopifnot(!is.null(outcomes))
      oi <- match(keys, paste(outcomes$patient_id, outcomes$documented_pregnancy_order, sep = ":"))
      qualifies <- !is.na(oi) & outcomes$recorded_live_birth[oi] %in% TRUE
      x[which(x$risk_code > 0 & !qualifies[idx]), risk_code := risk_code + 2L]
    }
  }
  if (specification == "postpartum_31_180") x[risk_code == 2L, risk_code := ifelse(pp_window == 2L, 2L, 3L)]
  if (specification == "postpartum_60day") x[risk_code == 2L, risk_code := 2L + (pp_window - 2L) %/% 2L]
  x <- x[included == TRUE]
  if (sccs) {
    retained <- x[, .(events = sum(events)), by = patient_id][events == 1L, patient_id]
    x <- x[patient_id %in% retained]
  }
  keys <- c("patient_id", "risk_code", "age_code", "calendar_code", if (art_duration) "art_duration_code")
  out <- x[, .(days = sum(days), events = sum(events)), by = keys]
  data.table::setorderv(out, keys)
  out
}

strict_tb_cohort <- function(cohort, evidence) {
  c <- data.table::copy(data.table::as.data.table(cohort))
  dates <- data.table::as.data.table(evidence)[source %in% c("treatment_start", "confirmed"),
    .(strict_date = min(evidence_date)), by = patient_id]
  c[, event_date := dates$strict_date[match(patient_id, dates$patient_id)]]
  c[!is.na(event_date) & (event_date < entry | event_date > administrative_exit), event_date := as.Date(NA)]
  c[, followup_end := data.table::fifelse(is.na(event_date), administrative_exit, event_date)]
  c
}

model_contrast <- function(model, numerator, denominator) {
  contrast <- setNames(rep(0, length(model$beta)), names(model$beta))
  contrast[numerator] <- 1; contrast[denominator] <- -1
  stopifnot(length(contrast) == length(model$beta))
  estimate <- sum(contrast * model$beta)
  se <- sqrt(as.numeric(t(contrast) %*% model$covariance %*% contrast))
  data.table::data.table(numerator = numerator, denominator = denominator, log_irr = estimate,
    standard_error = se, irr = exp(estimate), ci_lower = exp(estimate - 1.96 * se), ci_upper = exp(estimate + 1.96 * se),
    p_value = 2 * stats::pnorm(-abs(estimate / se)))
}

wald_contrast <- function(model, matrix, label) {
  effect <- as.numeric(matrix %*% model$beta)
  variance <- matrix %*% model$covariance %*% t(matrix)
  statistic <- as.numeric(t(effect) %*% solve(variance, effect))
  data.table::data.table(test = label, statistic = statistic, df = nrow(matrix),
    p_value = stats::pchisq(statistic, nrow(matrix), lower.tail = FALSE))
}

postpartum_tests <- function(model) {
  idx <- which(startsWith(names(model$beta), "postpartum_days_"))
  n <- length(idx); stopifnot(n >= 2L)
  selection <- diag(length(model$beta))[idx, , drop = FALSE]
  equality <- selection[-1L, , drop = FALSE] - matrix(selection[1, ], nrow = n - 1L, ncol = length(model$beta), byrow = TRUE)
  centered <- seq_len(n) - mean(seq_len(n)); weights <- centered / sum(centered^2)
  trend <- matrix(as.numeric(weights %*% selection), nrow = 1L)
  out <- data.table::rbindlist(list(wald_contrast(model, selection, "global_postpartum_vs_reference"),
    wald_contrast(model, equality, "heterogeneity_across_postpartum_windows"),
    wald_contrast(model, trend, "linear_trend_across_postpartum_windows")))
  estimate <- as.numeric(trend %*% model$beta)
  se <- sqrt(as.numeric(trend %*% model$covariance %*% t(trend)))
  window <- 180 / n
  out[test == "linear_trend_across_postpartum_windows", `:=`(
    window_days = window, log_rate_ratio_per_window_later = estimate, standard_error = se,
    rate_ratio_per_window_later = exp(estimate), ci_lower_per_window_later = exp(estimate - 1.96 * se),
    ci_upper_per_window_later = exp(estimate + 1.96 * se), rate_ratio_per_30_days_later = exp(estimate * 30 / window),
    ci_lower_per_30_days_later = exp((estimate - 1.96 * se) * 30 / window),
    ci_upper_per_30_days_later = exp((estimate + 1.96 * se) * 30 / window))]
  out
}
