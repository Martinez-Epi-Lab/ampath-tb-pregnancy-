#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "stage2_io.R", "eligibility_pregnancy.R", "person_time.R", "models.R", "secondary_analysis.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage2_secondary"); check_core_versions(project)
c <- stage2_input(ctx, "piece4", "cohort.rds"); p <- stage2_input(ctx, "piece4", "episodes.rds")
co <- stage2_input(ctx, "piece5", "cohort_intervals.rds"); sc <- stage2_input(ctx, "piece5", "sccs_intervals.rds")
b <- stage2_input(ctx, "piece6", "baseline_covariates.rds")
ev <- stage2_input(ctx, "piece4", "evidence.rds")
outcomes <- stage2_input(ctx, "stage2_outcomes", "outcomes.rds")
outcome_manifest <- jsonlite::read_json(file.path(project, "outputs/stage2_outcomes/run_manifest.json"))
cohort_manifest <- jsonlite::read_json(file.path(project, "outputs/piece4/run_manifest.json"))
stopifnot(identical(outcome_manifest$inputs$episodes$sha256, p$fingerprint$sha256),
  identical(outcome_manifest$inputs$visits$sha256, cohort_manifest$inputs$visits$sha256))
inputs <- list(cohort = c$fingerprint, episodes = p$fingerprint, cohort_intervals = co$fingerprint,
  sccs_intervals = sc$fingerprint, baseline = b$fingerprint, evidence = ev$fingerprint, outcomes = outcomes$fingerprint)
strict <- strict_tb_cohort(c$data, ev$data)
message("Constructing follow-up for the stricter TB definition.")
strict_time <- construct_person_time(strict, p$data)
outputs <- list(strict_cohort.rds = save_stage_rds(strict, ctx$work, "strict_cohort.rds"))
models <- descriptions <- contrasts <- timing_tests <- list()
specifications <- c("primary", "all_pregnancies", "order", "art_timing", "fully_observed", "live_birth", "postpartum_31_180", "postpartum_60day", "exclude_preconception", "utt_era", "cutoff_2024", "strict_tb")
for (spec in specifications) for (design in c("cohort", "sccs")) {
  id <- paste(spec, design, sep = "_")
  intervals <- if (spec == "strict_tb") strict_time[[paste0(design, "_intervals")]] else if (design == "cohort") co$data else sc$data
  rows <- secondary_rows(intervals, if (spec == "strict_tb") strict else c$data, p$data, spec, design == "sccs", outcomes$data)
  outputs[[paste0(id, "_rows.rds")]] <- save_stage_rds(rows, ctx$work, paste0(id, "_rows.rds"))
  labels <- secondary_labels(spec)
  variables <- c("age", if (spec == "utt_era") "calendar_utt" else "calendar")
  fn <- if (design == "cohort") fit_cohort_poisson else fit_sccs_poisson
  fit <- fn(rows, variables, labels = labels)
  fit$specification <- list(analysis = spec, design = design, variables = variables, labels = labels)
  models[[id]] <- fit
  description <- describe_person_time(rows, design, spec)
  description[, `:=`(model_id = id, feature = labels[risk_code + 1L])]
  descriptions[[id]] <- description
  if (spec %in% c("order", "art_timing")) {
    contrasts[[id]] <- rbindlist(list(model_contrast(fit, labels[4], labels[2]), model_contrast(fit, labels[5], labels[3])))
    contrasts[[id]][, model_id := id]
  }
  if (spec == "postpartum_60day") { timing_tests[[id]] <- postpartum_tests(fit); timing_tests[[id]][, model_id := id] }
  message(id, ": ", fit$diagnostics$analysis_n, " participants; ", fit$diagnostics$events, " events.")
}
art_rows <- secondary_rows(sc$data, c$data, p$data, "art_timing", sccs = TRUE, art_duration = TRUE)
id <- "art_timing_duration_sccs"
outputs[[paste0(id, "_rows.rds")]] <- save_stage_rds(art_rows, ctx$work, paste0(id, "_rows.rds"))
models[[id]] <- fit_sccs_poisson(art_rows, c("age", "calendar", "art_duration"), labels = secondary_labels("art_timing"))
models[[id]]$specification <- list(analysis = "art_timing", design = "sccs", variables = c("age", "calendar", "art_duration"),
  labels = secondary_labels("art_timing"))
descriptions[[id]] <- describe_person_time(art_rows, "sccs", "art_timing_duration")
descriptions[[id]][, `:=`(model_id = id, feature = secondary_labels("art_timing")[risk_code + 1L])]
contrasts[[id]] <- rbindlist(list(model_contrast(models[[id]], "pregnancy_conceived_after_art", "pregnancy_underway_at_art"),
  model_contrast(models[[id]], "postpartum_1_180_conceived_after_art", "postpartum_1_180_underway_at_art")))
contrasts[[id]][, model_id := id]
rich <- stage2_input(ctx, "piece5", "cohort_pooled_five.rds"); inputs$cohort_five <- rich$fingerprint
for (col in c("who_code", "cd4_code", "bmi_code", "clinic_code", "marital_code", "tpt_before_entry"))
  set(rich$data, j = col, value = b$data[[col]][match(rich$data$patient_id, b$data$patient_id)])
reference_variables <- c("age", "calendar", "who", "clinic", "marital")
influence_specs <- c(list(reference = reference_variables),
  setNames(lapply(reference_variables, function(v) setdiff(reference_variables, v)), paste0("without_", reference_variables)),
  setNames(lapply(c("tpt", "cd4", "bmi", "art_duration"), function(v) c(reference_variables, v)), paste0("plus_", c("tpt", "cd4", "bmi", "art_duration"))))
influence <- list()
for (nm in names(influence_specs)) {
  id <- paste0("influence_", nm)
  fit <- fit_cohort_poisson(rich$data, influence_specs[[nm]])
  fit$specification <- list(analysis = "covariate_influence", variables = influence_specs[[nm]])
  models[[id]] <- fit
  x <- copy(fit$estimates[feature %in% c("pregnancy", "postpartum_days_1_180")])
  ref <- models$influence_reference$estimates[match(x$feature, feature)]
  x[, `:=`(model_id = id, irr_change = irr - ref$irr, percent_irr_change = 100 * (irr / ref$irr - 1),
    percent_log_irr_change = 100 * (log_irr - ref$log_irr) / ref$log_irr)]
  influence[[id]] <- x
  message("Fitted ", id, ".")
}
core_models <- stage2_input(ctx, "piece6", "models.rds"); inputs$core_models <- core_models$fingerprint
for (design in c("cohort", "sccs")) {
  id <- paste0("postpartum_30day_", design)
  timing_tests[[id]] <- postpartum_tests(core_models$data[[paste0(design, "_30day")]])
  timing_tests[[id]][, model_id := id]
}
lag_contrasts <- rbindlist(lapply(c("cohort", "sccs"), function(design) {
  id <- paste0("postpartum_31_180_", design)
  result <- model_contrast(models[[id]], "postpartum_days_31_180", "postpartum_days_1_30")
  result[, model_id := id]; result
}))
fwrite(lag_contrasts, file.path(ctx$report, "postpartum_lag_contrasts.csv"))
fwrite(rbindlist(contrasts), file.path(ctx$report, "order_and_art_contrasts.csv"))
fwrite(rbindlist(timing_tests), file.path(ctx$report, "postpartum_interval_tests.csv"))
fwrite(rbindlist(influence), file.path(ctx$report, "covariate_influence.csv"))
outputs$models.rds <- save_stage_rds(models, ctx$work, "models.rds")
stage2_model_exports(ctx, models, descriptions)
summary <- list(models = length(models), cohort_women = nrow(c$data), strict_events = models$strict_tb_cohort$diagnostics$events,
  all_models_converged = all(vapply(models, function(m) m$diagnostics$converged, logical(1))))
finish_stage(ctx, inputs, outputs, list(specifications = specifications, influence_models = influence_specs), summary,
  c("R/secondary_analysis.R", "R/person_time.R", "R/models.R", "R/stage2_io.R"),
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, models_fitted = TRUE,
    first_two_restriction_applied = "Analysis-specific; all-pregnancy, order and ART-timing comparisons retain all orders"))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
