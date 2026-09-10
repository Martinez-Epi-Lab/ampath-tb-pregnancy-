#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "stage2_io.R", "pathway_mediation.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage3_mediation")
input <- stage2_input(ctx, "stage3_contact_cells", "cells.rds")
specifications <- data.table(model_id = c("primary", "expanded", "all_contacts", "known_routine", "lag_30day", "including_event_day"),
  mediator = c("routine_visit_dates", "routine_visit_dates", "all_visit_dates", "known_routine_visit_dates", "routine_visit_dates_30d_lag", "all_visit_dates_including_event_day"),
  covariates = c("age_calendar", "expanded", rep("age_calendar", 4)))
models <- effects <- coefficients <- diagnostics <- list()
for (i in seq_len(nrow(specifications))) {
  spec <- specifications[i]; id <- spec$model_id
  message("Fitting contact mediation: ", id, ".")
  model <- fit_contact_mediation(input$data, spec$mediator, spec$covariates)
  models[[id]] <- model
  effects[[id]] <- copy(model$effects); effects[[id]][, model_id := id]
  for (component in c("mediator", "total", "outcome")) {
    fit <- model[[component]]; se <- sqrt(diag(fit$covariance))
    coefficients[[paste(id, component)]] <- data.table(model_id = id, component = component, feature = names(fit$beta),
      coefficient = as.numeric(fit$beta), standard_error = se)
    diagnostics[[paste(id, component)]] <- as.data.table(c(list(model_id = id, component = component), fit$diagnostics))
    fwrite(data.table(feature = rownames(fit$covariance), as.data.table(fit$covariance)), file.path(ctx$report, paste0(id, "_", component, "_covariance.csv")))
  }
  fwrite(data.table(feature = rownames(model$joint_covariance), as.data.table(model$joint_covariance)), file.path(ctx$report, paste0(id, "_joint_covariance.csv")))
}
outputs <- list(models.rds = save_stage_rds(models, ctx$work, "models.rds"))
fwrite(rbindlist(effects), file.path(ctx$report, "mediation_effects.csv"))
fwrite(rbindlist(coefficients), file.path(ctx$report, "model_coefficients.csv"))
fwrite(rbindlist(diagnostics), file.path(ctx$report, "model_diagnostics.csv"))
summary <- list(specifications = nrow(specifications), component_fits = nrow(specifications) * 3L,
  exposure_result_rows = nrow(specifications) * 2L, patients = uniqueN(input$data$patient_id),
  events = sum(input$data$events), person_days = sum(input$data$days), all_models_converged = all(rbindlist(diagnostics)$converged))
finish_stage(ctx, list(cells = input$fingerprint), outputs, list(specifications = specifications,
  mediator_transformation = "log(1 + contacts per 100 observed days)", mediator_weight = "person-days",
  covariance_correction = "G/(G-1)", uncertainty = "Joint clustered multivariate delta method"), summary, "R/pathway_mediation.R",
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, first_two_restriction_applied = TRUE, models_fitted = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Mediation complete: six specifications and twelve exposure results.")
