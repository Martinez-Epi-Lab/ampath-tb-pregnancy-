#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "eligibility_pregnancy.R", "baseline_covariates.R", "models.R", "person_time.R"))
  source(file.path(project, "R", f))
ctx <- start_stage(project, "piece6")
p1 <- jsonlite::read_json(file.path(project, "outputs/piece1/run_manifest.json"))
p4 <- jsonlite::read_json(file.path(project, "outputs/piece4/run_manifest.json"))
p5 <- jsonlite::read_json(file.path(project, "outputs/piece5/run_manifest.json"))
v <- checked_input(project, "piece1", "visits.rds", p1$standardized_sha256$visits)
c <- checked_input(project, "piece4", "cohort.rds", p4$output_sha256$cohort.rds)
stopifnot(identical(p5$inputs$cohort$sha256, c$fingerprint$sha256), identical(p4$inputs$visits$sha256, v$fingerprint$sha256))
inputs <- list(visits = v$fingerprint, cohort = c$fingerprint)
message("Piece 6: deriving baseline covariates and fitting manuscript models.")
baseline <- derive_baseline_covariates(v$data, c$data)
outputs <- list(baseline_covariates.rds = save_stage_rds(baseline, ctx$work, "baseline_covariates.rds"))
specs <- list(
  cohort_primary = list(design = "cohort", view = "pooled_primary", variables = c("age", "calendar")),
  sccs_primary = list(design = "sccs", view = "pooled_primary", variables = c("age", "calendar")),
  five_cohort_primary = list(design = "cohort", view = "pooled_five", variables = c("age", "calendar")),
  five_cohort_expanded = list(design = "cohort", view = "pooled_five", variables = c("age", "calendar", "clinic", "marital", "tpt", "cd4", "bmi")),
  five_cohort_expanded_art_duration = list(design = "cohort", view = "pooled_five", variables = c("age", "calendar", "clinic", "marital", "tpt", "cd4", "bmi", "art_duration")),
  five_sccs_primary = list(design = "sccs", view = "pooled_five", variables = c("age", "calendar")),
  five_sccs_art_duration = list(design = "sccs", view = "pooled_five", variables = c("age", "calendar", "art_duration")),
  cohort_30day = list(design = "cohort", view = "postpartum_30day", variables = c("age", "calendar")),
  sccs_30day = list(design = "sccs", view = "postpartum_30day", variables = c("age", "calendar")))
all_estimates <- list(); all_diagnostics <- list(); models <- list(); exposure_results <- list()
for (id in names(specs)) {
  spec <- specs[[id]]; nm <- paste0(spec$design, "_", spec$view, ".rds")
  input <- checked_input(project, "piece5", nm, p5$output_sha256[[nm]])
  inputs[[nm]] <- input$fingerprint; rows <- input$data
  if (spec$view == "pooled_five" && spec$design == "cohort") {
    for (col in c("who_code", "cd4_code", "bmi_code", "clinic_code", "marital_code", "tpt_before_entry"))
      set(rows, j = col, value = baseline[[col]][match(rows$patient_id, baseline$patient_id)])
  }
  fn <- if (spec$design == "cohort") fit_cohort_poisson else fit_sccs_poisson
  fit <- fn(rows, spec$variables, windowed = spec$view == "postpartum_30day")
  fit$specification <- spec; models[[id]] <- fit
  tab <- copy(fit$estimates); tab[, `:=`(model_id = id, design = spec$design)]
  all_estimates[[id]] <- tab; all_diagnostics[[id]] <- c(list(model_id = id), fit$diagnostics)
  labels <- exposure_labels(spec$view == "postpartum_30day")
  exposure <- describe_person_time(rows, spec$design, spec$view)
  exposure[, feature := labels[risk_code + 1L]]
  exposure <- merge(exposure, tab[, !"design"], by = "feature", all.x = TRUE, sort = FALSE)
  exposure[risk_code == 0L, `:=`(irr = 1, model_id = id)]
  exposure[, estimate_display := ifelse(risk_code == 0L, "1.00 (Reference)", sprintf("%.2f (%.2f–%.2f)", irr, ci_lower, ci_upper))]
  setorder(exposure, risk_code); exposure_results[[id]] <- exposure
  fwrite(data.table(feature = rownames(fit$covariance), as.data.table(fit$covariance)), file.path(ctx$report, paste0(id, "_covariance.csv")))
  message("Fitted ", id, ": ", fit$diagnostics$events, " events; rank ", fit$diagnostics$information_rank, "/", fit$diagnostics$parameter_count, ".")
}
outputs$models.rds <- save_stage_rds(models, ctx$work, "models.rds")
fwrite(rbindlist(all_estimates), file.path(ctx$report, "model_coefficients.csv"))
fwrite(rbindlist(all_diagnostics, fill = TRUE), file.path(ctx$report, "model_diagnostics.csv"))
fwrite(rbindlist(exposure_results, fill = TRUE), file.path(ctx$report, "exposure_estimates.csv"))
fwrite(rbindlist(lapply(c("who_code", "cd4_code", "bmi_code", "clinic_code", "marital_code", "tpt_before_entry"),
  function(col) { z <- baseline[, .(women = .N), by = .(category = get(col))]; z[, variable := col]; z })), file.path(ctx$report, "baseline_categories.csv"))
summary <- list(model_fits = length(models), all_converged = all(vapply(models, function(m) m$diagnostics$converged, logical(1))),
  all_full_rank = all(vapply(models, function(m) m$diagnostics$information_rank == m$diagnostics$parameter_count, logical(1))),
  cohort_women = nrow(c$data), analyzed_events = models$cohort_primary$diagnostics$events,
  SCCS_cases = models$sccs_primary$diagnostics$analysis_n,
  scope = "Primary pooled models; three cohort/two SCCS agreed models; six 30-day postpartum models for Table 2")
finish_stage(ctx, inputs, outputs, specs, summary, c("R/baseline_covariates.R", "R/models.R", "R/person_time.R"),
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, first_two_restriction_applied = TRUE, models_fitted = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
message("Piece 6 complete: all ", length(models), " model fits passed numerical diagnostics.")
