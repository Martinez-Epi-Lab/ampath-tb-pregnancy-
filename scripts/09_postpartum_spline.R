#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "stage2_io.R", "eligibility_pregnancy.R", "models.R", "secondary_analysis.R", "postpartum_spline.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage2_spline"); check_core_versions(project)
p <- stage2_input(ctx, "piece4", "episodes.rds")
inputs <- list(episodes = p$fingerprint); outputs <- models <- curves <- tests <- list()
for (design in c("cohort", "sccs")) {
  input <- stage2_input(ctx, "piece5", paste0(design, "_intervals.rds")); inputs[[design]] <- input$fingerprint
  message("Constructing ", design, " postpartum spline rows.")
  rows <- spline_person_time(input$data, p$data, sccs = design == "sccs")
  outputs[[paste0(design, "_rows.rds")]] <- save_stage_rds(rows, ctx$work, paste0(design, "_rows.rds"))
  x <- spline_design_matrix(rows, intercept = design == "cohort")
  fn <- if (design == "cohort") fit_cohort_poisson else fit_sccs_poisson
  message("Fitting ", design, " spline: ", nrow(rows), " model rows.")
  fit <- fn(rows, x = x)
  fit$specification <- list(knots = c(10, 64, 117, 171), center = 90, design = design)
  models[[design]] <- fit
  curves[[design]] <- postpartum_spline_predictions(fit); curves[[design]][, design := design]
  tests[[design]] <- postpartum_spline_tests(fit); tests[[design]][, design := design]
  rm(rows, x, input); invisible(gc())
}
outputs$models.rds <- save_stage_rds(models, ctx$work, "models.rds")
stage2_model_exports(ctx, models, list())
fwrite(rbindlist(curves), file.path(ctx$report, "postpartum_spline_curve.csv"))
fwrite(rbindlist(tests), file.path(ctx$report, "postpartum_spline_tests.csv"))
summary <- list(models = length(models), all_converged = all(vapply(models, function(m) m$diagnostics$converged, logical(1))))
finish_stage(ctx, inputs, outputs, list(knots = c(10, 64, 117, 171), center = 90L, postpartum_days = 180L), summary,
  c("R/stage2_io.R", "R/models.R", "R/secondary_analysis.R", "R/postpartum_spline.R"),
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, first_two_restriction_applied = TRUE, models_fitted = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
