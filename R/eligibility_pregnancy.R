# Observation eligibility and pregnancy definitions.

study_rules <- function() list(
  data_start = as.Date("2010-01-01"), data_close = as.Date("2025-04-15"),
  run_in_days = 90L, entry_age_years = 15L, exit_age_years = 50L,
  max_pregnancy_elapsed_days = 315L
)

add_calendar_years <- function(x, years) {
  stopifnot(inherits(x, "Date"), length(years) == 1L, is.finite(years), years == floor(years))
  # Use calendar birthdays, with February 29 moved to February 28 in nonleap years.
  yr <- as.integer(format(x, "%Y")) + years
  md <- format(x, "%m-%d")
  target <- sprintf("%04d-%s", yr, md)
  target[is.na(x)] <- NA_character_
  out <- as.Date(target, format = "%Y-%m-%d")
  leap_clip <- !is.na(x) & md == "02-29" & is.na(out)
  out[leap_clip] <- as.Date(sprintf("%04d-02-28", yr[leap_clip]))
  if (any(!is.na(x) & is.na(out))) stop("Could not calculate calendar-year offset.")
  out
}

assert_unique_ids <- function(x, name) {
  if (anyNA(x$patient_id) || anyDuplicated(x$patient_id))
    stop(name, " requires unique, complete patient IDs.")
}

observation_windows <- function(demographics, visits, rules = study_rules()) {
  needed <- c("patient_id", "source_row", "male", "dob", "dod", "arvstartTHERAP")
  stopifnot(all(needed %in% names(demographics)), all(c("patient_id", "apptdate") %in% names(visits)))
  assert_unique_ids(demographics, "Demographics")
  for (nm in c("dob", "dod", "arvstartTHERAP")) stopifnot(inherits(demographics[[nm]], "Date"))
  stopifnot(inherits(visits$apptdate, "Date"))
  d <- data.table::copy(data.table::as.data.table(demographics)[, ..needed])
  data.table::setnames(d, "source_row", "demographics_source_row")
  v <- data.table::as.data.table(visits)
  bounds <- v[!is.na(patient_id) & !is.na(apptdate),
    .(first_visit = min(apptdate), last_visit = max(apptdate), valid_visit_records = .N), by = patient_id]
  vi <- match(d$patient_id, bounds$patient_id)
  data.table::set(d, j = "first_visit", value = bounds$first_visit[vi])
  data.table::set(d, j = "last_visit", value = bounds$last_visit[vi])
  data.table::set(d, j = "valid_visit_records", value = bounds$valid_visit_records[vi])
  data.table::set(d, j = "birthday_15", value = add_calendar_years(d$dob, rules$entry_age_years))
  data.table::set(d, j = "birthday_50", value = add_calendar_years(d$dob, rules$exit_age_years))
  reason <- rep("eligible_observation", nrow(d))
  mark <- function(condition, label) {
    idx <- which(reason == "eligible_observation" & !is.na(condition) & condition)
    reason[idx] <<- label
  }
  mark(is.na(d$first_visit) | is.na(d$last_visit), "no_visit_dates")
  mark(is.na(d$dob) | is.na(d$arvstartTHERAP), "missing_dob_or_art_start")
  mark(d$arvstartTHERAP <= rules$data_start, "art_not_after_2010_01_01")
  entry <- pmax(d$arvstartTHERAP, d$first_visit + rules$run_in_days, d$birthday_15, rules$data_start)
  exit <- pmin(d$last_visit, d$birthday_50, rules$data_close)
  has_death <- !is.na(d$dod)
  exit[has_death] <- pmin(exit[has_death], d$dod[has_death])
  not_calculated <- reason != "eligible_observation"
  entry[not_calculated] <- exit[not_calculated] <- as.Date(NA)
  mark(exit < entry, "exit_before_entry")
  data.table::set(d, j = "entry", value = entry)
  data.table::set(d, j = "exit", value = exit)
  data.table::set(d, j = "observation_status", value = reason)
  data.table::set(d, j = "female_coded", value = !is.na(d$male) & d$male == 0)
  # Record all entry and exit criteria that fall on the same date.
  data.table::set(d, j = "entry_at_art", value = !is.na(entry) & entry == d$arvstartTHERAP)
  data.table::set(d, j = "entry_at_run_in", value = !is.na(entry) & entry == d$first_visit + rules$run_in_days)
  data.table::set(d, j = "entry_at_15th_birthday", value = !is.na(entry) & entry == d$birthday_15)
  data.table::set(d, j = "exit_at_last_visit", value = !is.na(exit) & exit == d$last_visit)
  data.table::set(d, j = "exit_at_50th_birthday", value = !is.na(exit) & exit == d$birthday_50)
  data.table::set(d, j = "exit_at_data_close", value = !is.na(exit) & exit == rules$data_close)
  data.table::set(d, j = "exit_at_death", value = !is.na(exit) & !is.na(d$dod) & exit == d$dod)
  eligible <- d[observation_status == "eligible_observation"]
  data.table::set(eligible, j = "observation_days", value = as.integer(eligible$exit - eligible$entry) + 1L)
  stopifnot(all(eligible$observation_days >= 1L))
  list(audit = d, observation = eligible, female_observation = eligible[female_coded == TRUE], visit_bounds = bounds)
}

clean_pregnancy_episodes <- function(pregnancies, demographics, observation, rules = study_rules()) {
  needed <- c("patient_id", "source_row", "pregnancy", "intstart", "intend", "intdur")
  stopifnot(all(needed %in% names(pregnancies)), inherits(pregnancies$intstart, "Date"), inherits(pregnancies$intend, "Date"))
  assert_unique_ids(demographics, "Demographics")
  assert_unique_ids(observation, "Observation")
  p <- data.table::copy(data.table::as.data.table(pregnancies)[, ..needed])
  data.table::setnames(p, "source_row", "pregnancy_source_row")
  if (anyDuplicated(p$pregnancy_source_row) || anyNA(p$pregnancy_source_row)) stop("Pregnancy row provenance must be unique.")
  di <- match(p$patient_id, demographics$patient_id)
  oi <- match(p$patient_id, observation$patient_id)
  entry <- observation$entry[oi]; exit <- observation$exit[oi]
  original_start <- p$intstart; original_end <- p$intend
  duration <- as.integer(original_end - original_start)
  age <- as.numeric(original_start - demographics$dob[di]) / 365.25
  female <- !is.na(demographics$male[di]) & demographics$male[di] == 0
  status <- ifelse(!is.na(p$pregnancy) & p$pregnancy == 1, "candidate", "not_pregnancy_record")
  mark <- function(condition, label) {
    idx <- which(status == "candidate" & !is.na(condition) & condition)
    status[idx] <<- label
  }
  mark(is.na(di) | is.na(oi), "not_in_eligible_observation")
  mark(is.na(original_start) | is.na(original_end), "missing_start_or_end")
  mark(duration < 0, "negative_duration")
  mark(duration == 0, "zero_duration")
  mark(duration > rules$max_pregnancy_elapsed_days, "duration_over_45_weeks")
  mark(!female, "not_female_coded")
  mark(age < 15 | age >= 50, "pregnancy_start_age_outside_15_to_49")
  mark(original_end < entry, "ended_before_study_entry")
  mark(original_start > exit, "started_after_study_exit")
  data.table::set(p, j = "original_start", value = original_start)
  data.table::set(p, j = "original_end", value = original_end)
  data.table::set(p, j = "original_elapsed_days", value = duration)
  data.table::set(p, j = "age_at_original_start", value = age)
  data.table::set(p, j = "entry", value = entry)
  data.table::set(p, j = "exit", value = exit)
  candidate_rows <- which(status == "candidate")
  # Sort pregnancies by start, end, then source row for ties.
  idx <- candidate_rows[order(p$patient_id[candidate_rows], original_start[candidate_rows],
    original_end[candidate_rows], p$pregnancy_source_row[candidate_rows], method = "radix")]
  start <- original_start; end <- original_end
  overlap_adjusted <- rep(FALSE, nrow(p))
  shortened_by <- rep(NA_integer_, nrow(p))
  adjustments <- 0L; collapses <- 0L
  stack <- integer(length(idx)); top <- 0L; current_id <- NA_character_
  retained <- integer(length(idx)); retained_n <- 0L
  for (i in idx) {
    if (is.na(current_id) || p$patient_id[i] != current_id) {
      if (top > 0L) {
        retained[seq.int(retained_n + 1L, retained_n + top)] <- stack[seq_len(top)]
        retained_n <- retained_n + top
      }
      top <- 0L; current_id <- p$patient_id[i]
    }
    if (top > 0L && start[i] <= end[stack[top]]) {
      prior <- stack[top]
      adjustments <- adjustments + 1L
      end[prior] <- start[i] - 1L
      overlap_adjusted[prior] <- TRUE
      shortened_by[prior] <- p$pregnancy_source_row[i]
      if (end[prior] < start[prior]) {
        status[prior] <- "collapsed_after_overlap"
        top <- top - 1L; collapses <- collapses + 1L
      }
    }
    top <- top + 1L; stack[top] <- i
  }
  if (top > 0L) {
    retained[seq.int(retained_n + 1L, retained_n + top)] <- stack[seq_len(top)]
    retained_n <- retained_n + top
  }
  retained <- if (retained_n) retained[seq_len(retained_n)] else integer()
  analysis_start <- rep(as.Date(NA), nrow(p)); analysis_end <- analysis_start
  analysis_start[retained] <- pmax(start[retained], entry[retained])
  analysis_end[retained] <- pmin(end[retained], exit[retained])
  empty <- retained[analysis_end[retained] < analysis_start[retained]]
  status[retained] <- "retained"
  status[empty] <- "empty_after_observation_clipping"
  data.table::set(p, j = "start", value = start)
  data.table::set(p, j = "end", value = end)
  data.table::set(p, j = "analysis_start", value = analysis_start)
  data.table::set(p, j = "analysis_end", value = analysis_end)
  data.table::set(p, j = "left_truncated", value = original_start < entry)
  data.table::set(p, j = "right_truncated", value = original_end > exit)
  data.table::set(p, j = "overlap_adjusted", value = overlap_adjusted)
  data.table::set(p, j = "shortened_by_source_row", value = shortened_by)
  data.table::set(p, j = "pregnancy_status", value = status)
  clean <- p[pregnancy_status == "retained"]
  data.table::setorderv(clean, c("patient_id", "start", "end", "pregnancy_source_row"))
  clean[, documented_pregnancy_order := seq_len(.N), by = patient_id]
  clean[, first_two_documented := documented_pregnancy_order <= 2L]
  clean[, observed_pregnancy_days := as.integer(analysis_end - analysis_start) + 1L]
  if (nrow(clean)) {
    stopifnot(all(clean$observed_pregnancy_days >= 1L), all(clean$analysis_start >= clean$entry),
      all(clean$analysis_end <= clean$exit))
    overlap_check <- clean[, .(bad = any(start <= data.table::shift(end), na.rm = TRUE)), by = patient_id]
    if (any(overlap_check$bad)) stop("Cleaned pregnancy episodes still overlap.")
  }
  women <- observation[patient_id %in% clean$patient_id]
  if (nrow(women)) stopifnot(all(!is.na(women$male) & women$male == 0))
  list(audit = p, episodes = clean, pregnancy_eligible_observation = women,
    overlap_adjustments = adjustments, overlap_collapses = collapses,
    candidates_before_overlap = length(candidate_rows))
}
