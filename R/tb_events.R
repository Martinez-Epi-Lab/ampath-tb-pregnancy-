# TB evidence, event dates, and previous-TB exclusions.
tb_rules <- function() list(treatment_start_min = as.Date("2001-01-01"),
  treatment_start_max = as.Date("2025-04-15"), corroboration_days = 90L,
  prior_treatment = "Any nonmissing tbtxprior", prior_evidence = "Broad date strictly before entry")

extract_tb_evidence <- function(visits, rules = tb_rules()) {
  v <- data.table::as.data.table(visits)
  yes <- function(x) !is.na(x) & x == 1
  upper <- function(x) toupper(as.character(x))
  pc1 <- upper(v$tbpulmonaryPC); pc2 <- upper(v$tbextrapulmonaryPC)
  confirmed_pc <- pc1 %in% "C" | pc2 %in% "C"
  diagnosis <- yes(v$tbpulmdiag) | yes(v$tbextrapulmdiag) |
    yes(v$tbpulmonary) | yes(v$tbextrapulmonary) |
    pc1 %in% c("C", "P") | pc2 %in% c("C", "P")
  positive <- function(x) !is.na(x) & startsWith(upper(x), "POSITIVE")
  masks <- list(active = yes(v$activeTB), treatment = yes(v$tbtx),
    diagnosis = diagnosis, confirmed = confirmed_pc | positive(v$ZNstain) | positive(v$AFBculture),
    patient_type = !is.na(v$typetbpt),
    treatment_start = !is.na(v$tbtxstartdate) & v$tbtxstartdate >= rules$treatment_start_min &
      v$tbtxstartdate <= rules$treatment_start_max)
  rows <- lapply(names(masks), function(nm) {
    dates <- if (nm == "treatment_start") v$tbtxstartdate else v$apptdate
    idx <- which(masks[[nm]] & !is.na(v$patient_id) & !is.na(dates))
    data.table::data.table(patient_id = v$patient_id[idx], evidence_date = dates[idx],
      source = nm, visit_source_row = v$source_row[idx], appointment_date = v$apptdate[idx])
  })
  out <- data.table::rbindlist(rows)
  data.table::setorderv(out, c("patient_id", "evidence_date", "source", "visit_source_row"))
  out
}

pair_tb_dates <- function(first, second, gap = 90L) {
  # Evidence can occur in either order within the specified gap.
  # Use the earlier date for each qualifying pair.
  first <- sort(unique(as.integer(first))); second <- sort(unique(as.integer(second)))
  chunks <- lapply(first, function(a) {
    b <- second[abs(second - a) <= gap]
    if (!length(b)) return(NULL)
    data.table::data.table(candidate_date = pmin(a, b), evidence_date_1 = a, evidence_date_2 = b)
  })
  out <- data.table::rbindlist(chunks)
  if (!nrow(out)) out <- data.table::data.table(candidate_date = integer(), evidence_date_1 = integer(), evidence_date_2 = integer())
  for (nm in names(out)) data.table::set(out, j = nm, value = as.Date(out[[nm]], origin = "1970-01-01"))
  out
}

ascertain_tb <- function(evidence, rules = tb_rules()) {
  e <- unique(data.table::as.data.table(evidence)[, .(patient_id, source, evidence_date)])
  if (!nrow(e)) return(list(dates = data.table::data.table(patient_id = character(), broad_date = as.Date(character()),
    primary_date = as.Date(character()), confirmed_date = as.Date(character())),
    candidates = data.table::data.table(candidate_date = as.Date(character()), evidence_date_1 = as.Date(character()),
      evidence_date_2 = as.Date(character()), rule = character(), patient_id = character())))
  groups <- split(e, e$patient_id)
  candidates <- data.table::rbindlist(lapply(groups, function(g) {
    dates <- function(nm) g[source %in% nm, evidence_date]
    standalone <- g[source %in% c("treatment_start", "confirmed"),
      .(candidate_date = evidence_date, evidence_date_1 = evidence_date,
        evidence_date_2 = as.Date(NA), rule = source)]
    pairs <- list(diagnosis_active_or_treatment = list(dates("diagnosis"), dates(c("active", "treatment"))),
      active_treatment = list(dates("active"), dates("treatment")),
      patient_type_active_or_treatment = list(dates("patient_type"), dates(c("active", "treatment"))))
    linked <- lapply(names(pairs), function(nm) {
      z <- pair_tb_dates(pairs[[nm]][[1]], pairs[[nm]][[2]], rules$corroboration_days)
      z[, rule := rep(nm, .N)]; z
    })
    out <- data.table::rbindlist(c(list(standalone), linked), use.names = TRUE)
    out[, patient_id := rep(g$patient_id[1], .N)]
    out
  }), use.names = TRUE)
  if (!nrow(candidates)) candidates <- data.table::data.table(candidate_date = as.Date(character()),
    evidence_date_1 = as.Date(character()), evidence_date_2 = as.Date(character()), rule = character(), patient_id = character())
  data.table::setorderv(candidates, c("patient_id", "candidate_date", "rule", "evidence_date_1", "evidence_date_2"))
  dates <- e[, .(broad_date = min(evidence_date)), by = patient_id]
  first_primary <- candidates[, .(primary_date = min(candidate_date)), by = patient_id]
  first_confirmed <- e[source == "confirmed", .(confirmed_date = min(evidence_date)), by = patient_id]
  dates <- merge(merge(dates, first_primary, by = "patient_id", all.x = TRUE), first_confirmed, by = "patient_id", all.x = TRUE)
  list(dates = dates, candidates = candidates)
}

construct_tb_cohort <- function(observation, demographics, episodes, tb_dates) {
  audit <- data.table::copy(data.table::as.data.table(observation))
  stopifnot(!anyDuplicated(audit$patient_id), !anyDuplicated(tb_dates$patient_id))
  audit[, tbtxprior := demographics$tbtxprior[match(patient_id, demographics$patient_id)]]
  for (nm in c("broad_date", "primary_date", "confirmed_date"))
    data.table::set(audit, j = nm, value = tb_dates[[nm]][match(audit$patient_id, tb_dates$patient_id)])
  audit[, cohort_status := data.table::fcase(!is.na(tbtxprior), "tb_treatment_prior_to_enrollment",
    !is.na(broad_date) & broad_date < entry, "tb_evidence_before_study_entry", default = "included")]
  audit[, event_date := primary_date]
  audit[is.na(primary_date) | primary_date < entry | primary_date > exit | cohort_status != "included", event_date := as.Date(NA)]
  audit[, administrative_exit := exit]
  audit[, followup_end := data.table::fifelse(is.na(event_date), administrative_exit, event_date)]
  cohort <- audit[cohort_status == "included"]
  e <- data.table::copy(data.table::as.data.table(episodes)[patient_id %in% cohort$patient_id])
  list(audit = audit, cohort = cohort, episodes = e)
}
