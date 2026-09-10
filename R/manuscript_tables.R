# Central formatting for the manuscript tables.
format_number <- function(x, digits = 0L, signed = FALSE) {
  power <- 10^digits
  rounded <- sign(x) * floor(abs(x) * power + 0.5 + 1e-10) / power
  result <- formatC(rounded, format = "f", digits = digits, big.mark = ",")
  if (signed) result <- ifelse(rounded > 0, paste0("+", result), result)
  result[is.na(x)] <- ""
  result
}
format_interval <- function(value, lower, upper, digits = 2L, signed = FALSE, separator = "–")
  paste0(format_number(value, digits, signed), " (", format_number(lower, digits, signed), separator,
    format_number(upper, digits, signed), ")")
format_irr <- function(row) {
  stopifnot(nrow(row) == 1L)
  if ("risk_code" %in% names(row) && row$risk_code == 0L) "1.00 (Reference)" else format_interval(row$irr, row$ci_lower, row$ci_upper)
}
format_n_percent <- function(n, percent) paste0(format_number(n), " (", format_number(percent, 1), "%)")
format_iqr <- function(median, q1, q3, digits = 0L) paste0(format_number(median, digits), " (", format_number(q1, digits), "–", format_number(q3, digits), ")")
poisson_incidence <- function(events, days) {
  z <- stats::qnorm(0.975); y <- days / 365.25
  lower <- ifelse(events == 0L, 0, events * (1 - 1 / (9 * events) - z / (3 * sqrt(events)))^3)
  upper <- (events + 1) * (1 - 1 / (9 * (events + 1)) + z / (3 * sqrt(events + 1)))^3
  data.table::data.table(incidence = events / y * 100000, ci_lower = lower / y * 100000, ci_upper = upper / y * 100000)
}

build_manuscript_tables <- function(inputs) {
  core <- inputs$core; secondary <- inputs$secondary
  pick <- function(data, id, feature) { idx <- which(data$model_id == id & data$feature == feature); out <- data[idx]; if (nrow(out) != 1L) stop("Missing or duplicate result: ", id, " / ", feature); out }
  primary <- function(design, feature) pick(core, paste0(design, "_primary"), feature)
  period_label <- function(feature) ifelse(feature == "pregnancy", "Pregnancy", ifelse(feature == "other_observed_time", "Non-pregnant/non-postpartum time", gsub("([0-9]+) ([0-9]+)$", "\\1–\\2", gsub("_", " ", sub("^postpartum", "Postpartum", feature)))))
  tables <- list()
  add <- function(id, title, rows, notes, landscape = FALSE, groups = integer()) {
    rows <- as.data.frame(rows, check.names = FALSE, stringsAsFactors = FALSE)
    rows[] <- lapply(rows, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
    tables[[id]] <<- list(id = id, title = title, data = rows, notes = notes, landscape = landscape, groups = groups)
  }
  common <- "IRR, incidence rate ratio; SCCS, self-controlled case series; TB, tuberculosis. Pregnancy and postpartum are separate exposure categories; non-pregnant/non-postpartum time is the reference. Primary analyses include the first two documented pregnancies."
  add("Table_1", "Table 1  Characteristics of the cohort and primary SCCS population", inputs$table1,
    "Values are n (%) or median [IQR]. SCCS cases are nested within the cohort. Age and ART duration are measured at study entry. WHO stage, CD4 and BMI use the latest valid value in the 365 days before or on entry. Missing values are retained. Contributing pregnancies and follow-up use each design's observation rules; documented pregnancy order is not lifetime parity.",
    groups = which(inputs$table1[[2]] == "" & nzchar(inputs$table1[[1]])))
  rows <- list(); groups <- integer()
  for (design in c("cohort", "sccs")) {
    groups <- c(groups, length(rows) + 1L); rows[[length(rows) + 1L]] <- c(ifelse(design == "cohort", "Cohort", "SCCS"), "", "", "")
    selected_id <- paste0(design, "_30day"); x <- core[model_id == selected_id][order(risk_code)]
    for (i in seq_len(nrow(x))) rows[[length(rows) + 1L]] <- c(period_label(x$feature[i]), format_number(x$events[i]), format_number(x$person_years[i], 1), format_irr(x[i]))
  }
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Exposure interval", "TB events", "Person-years", "Adjusted IRR (95% CI)")
  add("Table_2", "Table 2  Tuberculosis risk during pregnancy and postpartum intervals", frame,
    paste(common, "Both models adjust for age and calendar period. Cohort follow-up stops at the first TB event; SCCS observation continues after the event. SCCS person-years are case observation time, not population incidence denominators."), groups = groups)
  specs <- c("primary", "all_pregnancies", "live_birth", "art_timing", "fully_observed", "postpartum_31_180", "exclude_preconception", "utt_era", "cutoff_2024")
  labels <- c("Primary analysis: first two documented pregnancies", "All documented pregnancies", "Pregnancies ending in a recorded live birth", "Pregnancies conceived after ART initiation", "Pregnancies fully observed after cohort entry", "Postpartum timing: days 31–180", "Reference excludes days 180–1 before pregnancy", "Observation from July 1, 2016", "Censoring on December 31, 2024")
  features <- list(live_birth = c("pregnancy_live_birth", "postpartum_live_birth"), art_timing = c("pregnancy_conceived_after_art", "postpartum_1_180_conceived_after_art"),
    fully_observed = c("pregnancy", "postpartum"), postpartum_31_180 = c("pregnancy", "postpartum_days_31_180"))
  out <- lapply(seq_along(specs), function(i) {
    spec <- specs[i]; f <- features[[spec]]; if (is.null(f)) f <- c("pregnancy", "postpartum_days_1_180")
    c(labels[i], unlist(lapply(c("cohort", "sccs"), function(d) vapply(f, function(feature) format_irr(pick(secondary, paste(spec, d, sep = "_"), feature)), character(1)))))
  })
  frame <- as.data.frame(do.call(rbind, out)); names(frame) <- c("Analysis", "Cohort pregnancy IRR (95% CI)", "Cohort postpartum IRR (95% CI)", "SCCS pregnancy IRR (95% CI)", "SCCS postpartum IRR (95% CI)")
  add("Table_3", "Table 3  Sensitivity analyses of pregnancy and postpartum tuberculosis risk", frame,
    paste(common, "All models adjust for age and calendar period. All-pregnancy and ART-timing comparisons retain all documented orders. Nonqualifying first-two periods remain nuisance categories in live-birth and fully observed analyses. Days 1–30 are modeled separately in the later-postpartum analysis. Preconception reference days are removed. The later-era restriction applies to person-time and allows earlier ART initiation."), TRUE)
  ids <- c("five_cohort_primary", "five_cohort_expanded", "five_cohort_expanded_art_duration", "five_sccs_primary", "five_sccs_art_duration")
  adjustment <- c("Age; calendar period", "Age; calendar period; clinic; partner status; recorded TPT; CD4; BMI", "Age; calendar period; clinic; partner status; recorded TPT; CD4; BMI; ART duration", "Age; calendar period", "Age; calendar period; ART duration")
  out <- lapply(seq_along(ids), function(i) c(c("Cohort model 1", "Cohort model 2", "Cohort model 3", "SCCS model 1", "SCCS model 2")[i], adjustment[i],
    vapply(c("pregnancy", "postpartum_days_1_180"), function(f) format_irr(pick(core, ids[i], f)), character(1))))
  frame <- as.data.frame(do.call(rbind, out)); names(frame) <- c("Model", "Adjustment variables", "Pregnancy IRR (95% CI)", "Postpartum IRR (95% CI)")
  add("Table_S1", "Table S1  Covariate adjustment in the cohort and SCCS models", frame, paste(common, "Recorded TPT refers to tuberculosis preventive treatment recorded on or before entry. Clinical covariates include explicit missing categories."), TRUE)
  inf <- inputs$influence
  infids <- paste0("influence_", c("reference", "without_age", "without_calendar", "without_who", "without_clinic", "without_marital", "plus_tpt", "plus_cd4", "plus_bmi", "plus_art_duration"))
  labels <- c("Reference model", "Remove age", "Remove calendar period", "Remove WHO stage", "Remove clinic", "Remove partner status", "Add recorded TPT", "Add baseline CD4", "Add baseline BMI", "Add ART duration")
  rows <- lapply(seq_along(infids), function(i) c(labels[i], unlist(lapply(c("pregnancy", "postpartum_days_1_180"), function(f) {
    z <- pick(inf, infids[i], f)
    change <- function(value, digits) if (i == 1L) "—" else format_number(value, digits, TRUE)
    c(format_number(z$irr, 2), change(z$irr_change, 4), if (i == 1L) "—" else paste0(change(z$percent_irr_change, 2), "%"), if (i == 1L) "—" else paste0(change(z$percent_log_irr_change, 2), "%"))
  }))))
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Model change", "Pregnancy IRR", "IRR change", "% IRR change", "% log IRR change", "Postpartum IRR", "IRR change", "% IRR change", "% log IRR change")
  add("Table_S2", "Table S2  Influence of covariate adjustment", frame,
    "The reference cohort model includes age, calendar period, WHO stage, clinic and partner status. Each row removes one variable or adds one variable to that model. All changes are calculated relative to the reference; signed differences and percentages are shown, including small changes.", TRUE)
  rows <- list()
  for (design in c("cohort", "sccs")) for (f in c("other_observed_time", "pregnancy", "postpartum_days_1_180")) {
    z <- primary(design, f); incidence <- "Not applicable"
    if (design == "cohort") { q <- poisson_incidence(z$events, z$person_days); incidence <- format_interval(q$incidence, q$ci_lower, q$ci_upper, 0) }
    rows[[length(rows) + 1L]] <- c(ifelse(design == "cohort", "Cohort", "SCCS"), period_label(f), format_number(z$contributors), format_number(z$events), format_number(z$person_years, 1), incidence, format_irr(z))
  }
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Design", "Exposure interval", "Contributing participants", "TB events", "Person-years", "Crude incidence per 100,000 person-years (95% CI)", "Adjusted IRR (95% CI)")
  add("Table_S3", "Table S3  Primary pooled exposure periods", frame, paste(common, "Crude incidence intervals use the original Wilson–Hilferty Poisson approximation. SCCS time is not a population incidence denominator."), TRUE)
  rows <- list()
  for (d in c("cohort", "sccs")) for (f in c("pregnancy", "postpartum_days_1_180")) rows[[length(rows) + 1L]] <- c(ifelse(d == "cohort", "Cohort", "SCCS"), period_label(f), format_irr(pick(secondary, paste0("all_pregnancies_", d), f)), format_irr(primary(d, f)))
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Design", "Period", "All pregnancies IRR (95% CI)", "First two pregnancies IRR (95% CI)")
  add("Table_S4", "Table S4  All pregnancies compared with the first two pregnancies", frame, paste(common, "Event counts reflect the pregnancy-order restriction used in each model."))
  rows <- lapply(c("cohort", "sccs"), function(d) c(ifelse(d == "cohort", "Cohort", "SCCS"), unlist(lapply(c("pregnancy_first_documented", "pregnancy_second_documented"), function(f) {
    z <- pick(secondary, paste0("order_", d), f); c(format_number(z$contributors), format_number(z$events), format_irr(z))
  }))))
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Design", "First pregnancy participants", "First pregnancy TB events", "First pregnancy IRR (95% CI)", "Second pregnancy participants", "Second pregnancy TB events", "Second pregnancy IRR (95% CI)")
  add("Table_S5", "Table S5  First and second documented pregnancies", frame, "These age- and calendar-adjusted models retain all documented pregnancy orders. First- and second-pregnancy estimates are displayed; third-or-later periods remain in the fitted model. Documented pregnancy order is not lifetime parity.", TRUE)
  phaseids <- c("preconception_days_180_to_1", "pregnancy", "postpartum_days_1_180"); phaselabels <- c("Preconception", "Pregnancy", "Postpartum")
  rows <- lapply(seq_along(phaseids), function(i) {
    x <- inputs$coverage[phase == phaseids[i]]
    c(phaselabels[i], format_number(x$eligible_episodes[1]), format_number(x$eligible_women[1]), vapply(c("cd4", "cd4_percentage", "viral_load"), function(f) { z <- x[candidate == f]; format_n_percent(z$episodes_measured, z$coverage_percent) }, character(1)))
  })
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Period", "Contributing episodes", "Contributing women", "CD4 count available", "CD4 percentage available", "Viral load available")
  biomarker_note <- "Preconception is days 180–1 before pregnancy and postpartum is days 1–180 after pregnancy. Windows are clipped to observed follow-up and postpartum stops before the next pregnancy. Analyses include the first two pregnancies and stop follow-up at the first TB event. The latest valid measurement in each phase is used."
  add("Table_S6", "Table S6  Availability of biomarker measurements", frame, biomarker_note, TRUE)
  pair <- inputs$biomarker_pairs
  comparison_label <- function(x) c(preconception_to_pregnancy = "Preconception to pregnancy", pregnancy_to_postpartum = "Pregnancy to postpartum", preconception_to_postpartum = "Preconception to postpartum")[x]
  x <- pair[candidate %in% c("cd4", "cd4_percentage")]
  rows <- lapply(seq_len(nrow(x)), function(i) {
    z <- x[i]; digits <- ifelse(z$candidate == "cd4", 0L, 1L)
    c(comparison_label(z$comparison), ifelse(z$candidate == "cd4", "CD4 count, cells/µL", "CD4 percentage, %"), paste0(format_number(z$paired_episodes), " (", format_number(z$paired_women), ")"),
      format_iqr(z$earlier_median, z$earlier_q1, z$earlier_q3, digits), format_iqr(z$later_median, z$later_q1, z$later_q3, digits), format_interval(z$change_estimate, z$change_ci_lower, z$change_ci_upper, digits, TRUE, " to "))
  })
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Comparison", "Measure", "Paired episodes (women)", "Earlier median (IQR)", "Later median (IQR)", "Median paired change (95% CI)")
  paired_note <- "Paired analyses require conception on or after entry and a measurement in both phases. Confidence intervals use 5,000 bootstrap samples of women, retaining all eligible paired episodes. CD4-count changes are cells/µL and CD4-percentage changes are percentage points."
  add("Table_S7A", "Table S7  Paired biomarker changes\nPanel A  CD4 count and CD4 percentage", frame, paired_note, TRUE)
  x <- pair[candidate == "viral_load"]
  rows <- lapply(seq_len(nrow(x)), function(i) { z <- x[i]; c(comparison_label(z$comparison), paste0(format_number(z$paired_episodes), " (", format_number(z$paired_women), ")"),
    paste0(format_number(z$earlier_n), "/", format_number(z$paired_episodes), " (", format_number(z$earlier_percent, 1), "%)"),
    paste0(format_number(z$later_n), "/", format_number(z$paired_episodes), " (", format_number(z$later_percent, 1), "%)"),
    paste0(format_number(z$changed_to_nonsuppressed_n), " / ", format_number(z$changed_to_suppressed_n)), format_interval(z$change_estimate, z$change_ci_lower, z$change_ci_upper, 2, TRUE, " to ")) })
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Comparison", "Paired episodes (women)", "Earlier nonsuppression n/N (%)", "Later nonsuppression n/N (%)", "Became nonsuppressed / suppressed", "Paired prevalence difference (95% CI)")
  add("Table_S7B", "Table S7  Paired biomarker changes\nPanel B  Viral load nonsuppression", frame, "Viral nonsuppression is at least 1,000 copies/mL. Differences are percentage points. Bootstrap intervals resample women and retain paired episodes.", TRUE)
  x <- inputs$contact_rates
  rows <- lapply(seq_len(nrow(x)), function(i) { z <- x[i]; c(ifelse(z$exposure == "unexposed", "Not pregnant or postpartum", period_label(z$exposure)), format_number(z$patients), format_number(z$person_days), format_number(z$routine_visit_dates), format_interval(z$routine_visit_dates_per_100_days, z$routine_visit_rate_ci_lower, z$routine_visit_rate_ci_upper)) })
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Period", "Contributing women", "Observed person-days", "Routine non-TB contact dates", "Contacts per 100 days (95% CI)")
  contact_note <- "Multiple records on the same date count as one contact. Routine contacts exclude dates with TB encounter codes or flags, chest-radiography entries or bacteriological testing, and exclude the TB event day. Confidence intervals use 5,000 woman-clustered bootstrap samples."
  add("Table_S8", "Table S8  Routine healthcare-contact rates", frame, contact_note, TRUE)
  x <- inputs$contact_pairs
  rows <- lapply(seq_len(nrow(x)), function(i) { z <- x[i]; c(comparison_label(z$comparison), paste0(format_number(z$paired_episodes), " (", format_number(z$paired_women), ")"), format_iqr(z$earlier_median, z$earlier_q1, z$earlier_q3, 2), format_iqr(z$later_median, z$later_q1, z$later_q3, 2), format_interval(z$median_paired_change, z$change_ci_lower, z$change_ci_upper, 2, TRUE, " to ")) })
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Comparison", "Paired episodes (women)", "Earlier median rate (IQR)", "Later median rate (IQR)", "Median paired change (95% CI)")
  add("Table_S9", "Table S9  Paired changes in routine healthcare-contact rates", frame, paste(contact_note, "Rates and changes are contacts per 100 observed days. Paired analyses require conception on or after entry and observation in both phases."), TRUE)
  med <- inputs$mediation; x <- med[model_id %in% c("primary", "expanded")]
  rows <- lapply(seq_len(nrow(x)), function(i) { z <- x[i]; c(ifelse(z$model_id == "primary", "Age and calendar period", "Expanded adjustment"), period_label(z$exposure), vapply(c("total_effect", "direct_effect", "indirect_effect", "decomposed_total"), function(prefix) format_interval(z[[paste0(prefix, "_rr")]], z[[paste0(prefix, "_ci_lower")]], z[[paste0(prefix, "_ci_upper")]]), character(1))) })
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Adjustment", "Exposure", "Total-effect IRR (95% CI)", "Natural direct RR (95% CI)", "Natural indirect RR (95% CI)", "Decomposed total RR (95% CI)")
  mednote <- "The mediator is log(1 + routine contacts per 100 observed days). Models include exposure–mediator interactions. Effects are standardized over observed person-time, with joint woman-clustered covariance and delta-method intervals. Expanded adjustment adds ART duration, clinic, partner status, recorded TPT, CD4, BMI and pre-entry contact intensity. The separately fitted total effect differs from the decomposed total."
  add("Table_S10A", "Table S10  Healthcare-contact mediation\nPanel A  Primary and expanded adjustment", frame, mednote, TRUE)
  mids <- c("primary", "all_contacts", "known_routine", "lag_30day", "including_event_day")
  labels <- c("Routine non-TB contacts; event day excluded", "All contacts; event day excluded", "Routine contacts with a known encounter code", "Routine non-TB contacts at least 30 days before TB", "All contacts including the TB event day")
  rows <- lapply(seq_along(mids), function(i) c(labels[i], vapply(c("pregnancy", "postpartum_days_1_180"), function(f) { z <- med[model_id == mids[i] & exposure == f]; format_interval(z$indirect_effect_rr, z$indirect_effect_ci_lower, z$indirect_effect_ci_upper) }, character(1))))
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Contact definition", "Pregnancy natural indirect RR (95% CI)", "Postpartum natural indirect RR (95% CI)")
  add("Table_S10B", "Table S10  Healthcare-contact mediation\nPanel B  Contact-definition sensitivities", frame, "All sensitivity models adjust for age and calendar period. Routine contacts exclude dates containing TB codes, indicators or testing. The known-code analysis also requires an encounter code. Lagged contacts occur on or before TB day −30.", TRUE)
  rows <- list()
  for (d in c("cohort", "sccs")) for (f in c("pregnancy", "postpartum")) {
    z <- pick(secondary, paste0("live_birth_", d), paste0(f, "_live_birth")); basef <- ifelse(f == "pregnancy", f, "postpartum_days_1_180")
    rows[[length(rows) + 1L]] <- c(ifelse(d == "cohort", "Cohort", "SCCS"), ifelse(f == "pregnancy", "Pregnancy", "Postpartum"), format_irr(primary(d, basef)), format_irr(z), format_number(z$events))
  }
  frame <- as.data.frame(do.call(rbind, rows)); names(frame) <- c("Design", "Period", "Primary IRR (95% CI)", "Recorded live-birth IRR (95% CI)", "Exposed TB events")
  q <- inputs$outcome_summary
  add("Table_S11", "Table S11  Recorded-live-birth sensitivity analysis", frame, paste0(common, " Outcome documentation was present for ", format_number(q$primary_recorded_episodes), " of ", format_number(q$first_two_episodes), " episodes (", format_number(q$primary_recorded_episodes / q$first_two_episodes * 100, 1), "%); ", format_number(q$recorded_live_birth_episodes), " were classified as live births. Evidence was linked from 30 days before through 90 days after the original pregnancy end. Other first-two pregnancy-related time remained in separate nuisance categories."), TRUE)
  tables
}
