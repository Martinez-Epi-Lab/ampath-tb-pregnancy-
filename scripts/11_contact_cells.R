#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "stage2_io.R", "pathway_contacts.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage3_contact_cells"); check_core_versions(project)
v <- stage2_input(ctx, "piece1", "visits.rds")
c <- stage2_input(ctx, "piece4", "cohort.rds")
x <- stage2_input(ctx, "piece5", "cohort_intervals.rds")
b <- stage2_input(ctx, "piece6", "baseline_covariates.rds")
message("Classifying distinct healthcare-contact dates.")
day <- build_contact_dates(v$data, c$data)
message("Counting contacts within exposure, age, calendar and ART-duration cells.")
cells <- build_contact_cells(x$data, day, c$data, b$data)
outputs <- list(contact_dates.rds = save_stage_rds(day, ctx$work, "contact_dates.rds"),
  cells.rds = save_stage_rds(cells, ctx$work, "cells.rds"))
phase <- contact_phase_summary(cells)
fwrite(phase, file.path(ctx$report, "contact_phase_summary.csv"))
summary <- list(cohort_women = uniqueN(cells$patient_id), cells = nrow(cells), person_days = sum(cells$days),
  events = sum(cells$events), distinct_source_contact_dates = nrow(day), routine_contact_dates = sum(cells$routine_visit_dates))
finish_stage(ctx, list(visits = v$fingerprint, cohort = c$fingerprint, intervals = x$fingerprint, baseline = b$fingerprint),
  outputs, list(contact_definitions = contact_count_columns(), event_day_excluded_from_primary = TRUE,
    baseline_contact_lookback_days = 180L), summary, "R/pathway_contacts.R",
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, first_two_restriction_applied = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Contact cells complete: ", nrow(cells), " cells and ", summary$routine_contact_dates, " routine contact dates.")
