# Cohort and SCCS person-time. Interval endpoints are inclusive.
person_time_rules <- function() list(postpartum_days = 180L, postpartum_window_days = 30L,
  maximum_documented_pregnancies = 2L, days_per_year = 365.25,
  endpoints = "Inclusive", cohort_censoring = "First TB event, including excluded exposure time",
  SCCS = "Full observation for cases with their event in retained time")

person_intervals <- function(record, pregnancies) {
  r <- record; e <- pregnancies
  entry <- as.integer(r$entry); exit <- as.integer(r$administrative_exit)
  event <- as.integer(r$event_date)
  starts <- as.integer(e$analysis_start); ends <- as.integer(e$analysis_end)
  original_starts <- as.integer(e$start); adjusted_ends <- as.integer(e$end)
  n <- nrow(e)
  next_start <- c(original_starts[-1], Inf)
  pp_start <- adjusted_ends + 1L
  pp_end <- pmin(adjusted_ends + 180L, next_start - 1L, exit)
  birthdays <- vapply(c(20L, 30L, 40L), function(y) as.integer(add_calendar_years(r$dob, y)), integer(1))
  calendars <- as.integer(as.Date(c("2015-01-01", "2020-01-01")))
  arts <- as.integer(r$arvstartTHERAP) + c(180L, 365L, 730L, 1825L)
  pp_breaks <- unlist(lapply(adjusted_ends, function(d) d + c(1L, 31L, 61L, 91L, 121L, 151L, 181L)))
  breaks <- sort(unique(c(entry, exit + 1L, starts, ends + 1L, pp_start, pp_end + 1L,
    pp_breaks, birthdays, calendars, arts, event, event + 1L)))
  breaks <- breaks[is.finite(breaks) & breaks >= entry & breaks <= exit + 1L]
  left <- head(breaks, -1L); right <- tail(breaks, -1L) - 1L
  risk <- order <- pp <- integer(length(left))
  episode_order <- integer(length(left))
  for (i in seq_len(n)) {
    idx <- which(left >= pp_start[i] & left <= pp_end[i])
    risk[idx] <- 2L; order[idx] <- if (i <= 2L) 2L * i else 6L
    pp[idx] <- 2L + (left[idx] - pp_start[i]) %/% 30L
    episode_order[idx] <- i
  }
  for (i in seq_len(n)) {
    idx <- which(left >= starts[i] & left <= ends[i])
    risk[idx] <- 1L; order[idx] <- if (i <= 2L) 2L * i - 1L else 5L
    pp[idx] <- 1L; episode_order[idx] <- i
  }
  age <- ifelse(left < birthdays[1], 1L, ifelse(left < birthdays[2], 0L, ifelse(left < birthdays[3], 2L, 3L)))
  calendar <- findInterval(left, calendars)
  art <- findInterval(left, arts)
  data.table::data.table(patient_id = rep(r$patient_id, length(left)),
    interval_start = as.Date(left, origin = "1970-01-01"), interval_end = as.Date(right, origin = "1970-01-01"),
    risk_code = risk, order_code = order, pp_window = pp, documented_pregnancy_order = episode_order,
    age_code = age, calendar_code = calendar, art_duration_code = art,
    days = as.integer(right - left + 1L), events = as.integer(!is.na(event) & left <= event & right >= event),
    included = !(order %in% c(5L, 6L)))
}

construct_person_time <- function(cohort, episodes) {
  c <- data.table::copy(data.table::as.data.table(cohort))
  p <- data.table::copy(data.table::as.data.table(episodes))
  data.table::setorderv(c, "patient_id"); data.table::setorderv(p, c("patient_id", "documented_pregnancy_order"))
  grouped <- split(p, p$patient_id)
  all <- data.table::rbindlist(lapply(seq_len(nrow(c)), function(i) person_intervals(c[i], grouped[[c$patient_id[i]]])))
  event_audit <- all[events == 1L, .(patient_id, event_date = interval_start, risk_code, order_code, pp_window,
    included_in_primary_analysis = included)]
  cohort_rows <- all[interval_end <= c$followup_end[match(patient_id, c$patient_id)]]
  sccs_rows <- all[patient_id %in% c$patient_id[!is.na(c$event_date)]]
  retained_cases <- event_audit[included_in_primary_analysis == TRUE, patient_id]
  sccs_rows[, primary_case := patient_id %in% retained_cases]
  # Check that each observed day is counted once.
  conservation <- all[, .(actual_days = sum(days)), by = patient_id]
  stopifnot(all(conservation$actual_days == as.integer(c$administrative_exit[match(conservation$patient_id, c$patient_id)] -
    c$entry[match(conservation$patient_id, c$patient_id)]) + 1L), all(all$days > 0L))
  list(cohort_intervals = cohort_rows, sccs_intervals = sccs_rows, event_audit = event_audit)
}

analysis_rows <- function(intervals, view = c("pooled_primary", "pooled_five", "postpartum_30day"), sccs = FALSE) {
  view <- match.arg(view)
  x <- data.table::copy(data.table::as.data.table(intervals)[included == TRUE])
  if (sccs) x <- x[primary_case == TRUE]
  keys <- c("patient_id", "risk_code", "age_code", "calendar_code")
  if (view == "pooled_five") keys <- c(keys, "order_code", "art_duration_code")
  if (view == "postpartum_30day") x[, risk_code := pp_window]
  out <- x[, .(days = sum(days), events = sum(events)), by = keys]
  # Keep pregnancy-order and ART-duration cells separate in the five-model set.
  # The variance correction depends on the number of model rows.
  data.table::setorderv(out, if (view == "pooled_five") c("patient_id", "order_code", "age_code", "calendar_code", "art_duration_code") else keys)
  out
}

describe_person_time <- function(rows, design, view) {
  out <- rows[, .(events = sum(events), person_days = sum(days), contributors = data.table::uniqueN(patient_id)), by = risk_code]
  out[, `:=`(person_years = person_days / 365.25, design = design, view = view)]
  if (design == "cohort") out[, incidence_per_100000_person_years := events / person_years * 100000]
  data.table::setorder(out, risk_code); out
}
