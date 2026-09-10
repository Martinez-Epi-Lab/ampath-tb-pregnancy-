# Distinct contact dates and contact counts within the analysis cells.
contact_count_columns <- function() c("all_visit_dates", "routine_visit_dates", "known_routine_visit_dates",
  "routine_visit_dates_30d_lag", "all_visit_dates_including_event_day", "routine_visit_dates_including_event_day")

build_contact_dates <- function(visits, cohort) {
  flags <- c("activeTB", "tbtx", "tbpulmdiag", "tbextrapulmdiag", "tbpulmonary", "tbextrapulmonary",
    "tbpulmonaryPC", "tbextrapulmonaryPC", "lymphnodeTB", "lymphnodeTBDiag")
  fields <- c("patient_id", "apptdate", "encounter_type", "combined_enc_types", "cxr", "ZNstain", "AFBculture", flags)
  v <- data.table::copy(data.table::as.data.table(visits)[patient_id %in% cohort$patient_id & !is.na(apptdate), ..fields])
  combined <- gsub(" ", "", v$combined_enc_types, fixed = TRUE)
  coded <- v$encounter_type %in% c(27, 95) | (!is.na(combined) & grepl("(^|,)(27|95)(,|$)", combined))
  for (col in flags) coded <- coded | (!is.na(v[[col]]) & v[[col]] != 0)
  for (col in c("cxr", "ZNstain", "AFBculture")) coded <- coded | (!is.na(v[[col]]) & trimws(as.character(v[[col]])) != "")
  v[, `:=`(tb_coded_contact = coded, encounter_code_available = !is.na(encounter_type) | !is.na(combined_enc_types))]
  day <- v[, .(tb_coded_contact = any(tb_coded_contact), encounter_code_available = any(encounter_code_available)),
    by = .(patient_id, date = apptdate)]
  day[, `:=`(routine_non_tb_contact = !tb_coded_contact,
    known_routine_non_tb_contact = encounter_code_available & !tb_coded_contact)]
  data.table::setorder(day, patient_id, date)
  day
}

pre_entry_contact_rate <- function(day, cohort) {
  windows <- cohort[, .(patient_id, start = entry - 180L, end = entry - 1L)]
  counts <- day[windows, on = .(patient_id, date >= start, date <= end), .(contacts = .N), by = .EACHI, nomatch = 0L]
  result <- data.table::data.table(patient_id = cohort$patient_id)
  result[, contacts := counts$contacts[match(patient_id, counts$patient_id)]]
  result[is.na(contacts), contacts := 0L]
  result[, pre_entry_visit_rate_log := log1p(contacts / 180 * 100)]
  result
}

build_contact_cells <- function(intervals, day, cohort, baseline) {
  x <- data.table::copy(data.table::as.data.table(intervals)[included == TRUE])
  keys <- c("patient_id", "risk_code", "age_code", "calendar_code", "art_duration_code")
  cells <- x[, .(days = sum(days), events = sum(events)), by = keys]
  hits <- x[day, on = .(patient_id, interval_start <= date, interval_end >= date), nomatch = 0L,
    .(patient_id = x.patient_id, date = i.date, risk_code, age_code, calendar_code, art_duration_code,
      routine = i.routine_non_tb_contact, known = i.known_routine_non_tb_contact)]
  stopifnot(!anyDuplicated(hits[, .(patient_id, date)]))
  hits[, event_date := cohort$event_date[match(patient_id, cohort$patient_id)]]
  hits[, before_event := is.na(event_date) | date < event_date]
  counts <- hits[, .(all_visit_dates = sum(before_event), routine_visit_dates = sum(before_event & routine),
    known_routine_visit_dates = sum(before_event & known),
    routine_visit_dates_30d_lag = sum(before_event & routine & (is.na(event_date) | date <= event_date - 30L)),
    all_visit_dates_including_event_day = .N, routine_visit_dates_including_event_day = sum(routine)), by = keys]
  cells <- merge(cells, counts, by = keys, all.x = TRUE, sort = FALSE)
  for (col in contact_count_columns()) data.table::set(cells, which(is.na(cells[[col]])), col, 0L)
  for (col in c("cd4_code", "bmi_code", "clinic_code", "marital_code", "tpt_before_entry"))
    data.table::set(cells, j = col, value = baseline[[col]][match(cells$patient_id, baseline$patient_id)])
  pre <- pre_entry_contact_rate(day, cohort)
  cells[, pre_entry_visit_rate_log := pre$pre_entry_visit_rate_log[match(patient_id, pre$patient_id)]]
  data.table::setorderv(cells, keys)
  stopifnot(!anyNA(cells), sum(cells$days) == sum(x$days), sum(cells$events) == sum(x$events))
  cells
}

contact_phase_summary <- function(cells) {
  result <- cells[, c(list(patients = data.table::uniqueN(patient_id), person_days = sum(days), tb_events = sum(events)),
    lapply(.SD, sum)), by = .(risk_code), .SDcols = contact_count_columns()]
  result[, `:=`(exposure = c("unexposed", "pregnancy", "postpartum_days_1_180")[risk_code + 1L], person_years = person_days / 365.25)]
  result[, tb_incidence_per_100000_person_years := tb_events / person_years * 100000]
  for (col in contact_count_columns()) data.table::set(result, j = paste0(col, "_per_100_days"), value = result[[col]] / result$person_days * 100)
  data.table::setorder(result, risk_code)
  result
}

contact_episode_counts <- function(phases, day, cohort) {
  dates <- data.table::copy(day[routine_non_tb_contact == TRUE])
  dates[, event_date := cohort$event_date[match(patient_id, cohort$patient_id)]]
  dates <- dates[is.na(event_date) | date < event_date]
  windows <- data.table::copy(phases)
  windows[, phase_row := seq_len(.N)]
  # The interval key is unique within each woman's documented episode and phase.
  hit_rows <- dates[windows, on = .(patient_id, date >= start, date <= end), nomatch = 0L, .(phase_row = i.phase_row)]
  counts <- hit_rows[, .(routine_contact_dates = .N), by = phase_row]
  windows[, routine_contact_dates := counts$routine_contact_dates[match(phase_row, counts$phase_row)]]
  windows[is.na(routine_contact_dates), routine_contact_dates := 0L]
  windows[, `:=`(observed_days = as.integer(end - start) + 1L, phase_row = NULL)]
  windows[, contacts_per_100_days := routine_contact_dates / observed_days * 100]
  windows
}

contact_comparisons <- function() list(preconception_to_pregnancy = c("preconception_days_180_to_1", "pregnancy"),
  pregnancy_to_postpartum = c("pregnancy", "postpartum_days_1_180"),
  preconception_to_postpartum = c("preconception_days_180_to_1", "postpartum_days_1_180"))

contact_pair_frames <- function(phases) {
  lapply(contact_comparisons(), function(periods) {
    a <- phases[conceived_after_entry == TRUE & phase == periods[1],
      .(patient_id, pregnancy_index, episode_number, earlier_rate = contacts_per_100_days)]
    b <- phases[conceived_after_entry == TRUE & phase == periods[2],
      .(patient_id, pregnancy_index, later_rate = contacts_per_100_days)]
    z <- merge(a, b, by = c("patient_id", "pregnancy_index"))
    z[, change := later_rate - earlier_rate]
    data.table::setorder(z, episode_number)
    z
  })
}

bootstrap_contact_change <- function(pairs, rng, replicates = 5000L, chunk_size = 250L) {
  women <- unique(pairs$patient_id); women <- women[order(as.numeric(women))]
  code <- match(pairs$patient_id, women)
  ord <- order(pairs$change, method = "radix"); values <- pairs$change[ord]; code <- code[ord]
  estimates <- numeric(replicates)
  for (start in seq.int(1L, replicates, by = chunk_size)) {
    size <- min(chunk_size, replicates - start + 1L)
    weights <- manuscript_rng_counts(rng, length(women), size, "multinomial")
    for (j in seq_len(size)) {
      w <- weights[code, j]; cumulative <- cumsum(w); target <- (sum(w) - 1) / 2
      lo <- findInterval(floor(target), cumulative) + 1L
      hi <- findInterval(ceiling(target), cumulative) + 1L
      estimates[start + j - 1L] <- (values[lo] + values[hi]) / 2
    }
  }
  estimates
}

bootstrap_contact_rates <- function(cells, rng, replicates = 5000L, chunk_size = 250L) {
  women <- unique(cells$patient_id); women <- women[order(as.numeric(women))]
  pp <- cells[, .(days = sum(days), contacts = sum(routine_visit_dates)), by = .(patient_id, risk_code)]
  totals <- matrix(0, nrow = length(women), ncol = 6L)
  for (risk in 0:2) {
    z <- pp[risk_code == risk]; idx <- match(z$patient_id, women)
    totals[idx, 2L * risk + 1L] <- z$days; totals[idx, 2L * risk + 2L] <- z$contacts
  }
  estimates <- matrix(NA_real_, nrow = replicates, ncol = 3L)
  for (start in seq.int(1L, replicates, by = chunk_size)) {
    size <- min(chunk_size, replicates - start + 1L)
    weights <- manuscript_rng_counts(rng, length(women), size, "multinomial")
    sums <- crossprod(weights, totals)
    estimates[start + seq_len(size) - 1L, ] <- sums[, c(2, 4, 6)] / sums[, c(1, 3, 5)] * 100
  }
  colnames(estimates) <- c("unexposed", "pregnancy", "postpartum_days_1_180")
  estimates
}

summarize_contact_pairs <- function(pairs, draws, comparison) {
  q <- function(x) unname(stats::quantile(x, c(0.5, 0.25, 0.75), type = 7))
  before <- q(pairs$earlier_rate); after <- q(pairs$later_rate)
  ci <- unname(stats::quantile(draws, c(0.025, 0.975), type = 7))
  data.table::data.table(comparison = comparison, paired_episodes = nrow(pairs), paired_women = data.table::uniqueN(pairs$patient_id),
    earlier_median = before[1], earlier_q1 = before[2], earlier_q3 = before[3], later_median = after[1], later_q1 = after[2], later_q3 = after[3],
    median_paired_change = stats::median(pairs$change), change_ci_lower = ci[1], change_ci_upper = ci[2])
}
