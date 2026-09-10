#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
source(file.path(project, "R", "stage_io.R"))
source(file.path(project, "R", "eligibility_pregnancy.R"))
ctx <- start_stage(project, "piece3")
p1 <- jsonlite::read_json(file.path(project, "outputs", "piece1", "run_manifest.json"))
p2 <- jsonlite::read_json(file.path(project, "outputs", "piece2", "run_manifest.json"))
message("Piece 3: checking inputs and reconstructing pregnancy episodes.")
d <- checked_input(project, "piece1", "demographics.rds", p1$standardized_sha256$demographics)
p <- checked_input(project, "piece1", "pregnancies.rds", p1$standardized_sha256$pregnancies)
o <- checked_input(project, "piece2", "observation.rds", p2$output_sha256$observation.rds)
if (!identical(p2$inputs$demographics$sha256, d$fingerprint$sha256) ||
  !identical(p2$inputs$visits$sha256, p1$standardized_sha256$visits))
  stop("Piece 2 belongs to a different Piece 1 dataset version; rerun Piece 2.")
rules <- study_rules()
result <- clean_pregnancy_episodes(p$data, d$data, o$data, rules)
inputs <- list(demographics = d$fingerprint, pregnancies = p$fingerprint, observation = o$fingerprint)
outputs <- list()
for (nm in c("audit", "episodes", "pregnancy_eligible_observation"))
  outputs[[paste0(nm, ".rds")]] <- save_stage_rds(result[[nm]], ctx$work, paste0(nm, ".rds"))
flow <- result$audit[, .(interval_rows = .N, unique_people = uniqueN(patient_id)), by = pregnancy_status]
setorderv(flow, "pregnancy_status")
fwrite(flow, file.path(ctx$report, "pregnancy_cleaning_flow.csv"))
e <- result$episodes
order_counts <- e[, .(episodes = .N, people = uniqueN(patient_id), observed_pregnancy_days = sum(observed_pregnancy_days)), by = documented_pregnancy_order]
setorder(order_counts, documented_pregnancy_order)
fwrite(order_counts, file.path(ctx$report, "documented_pregnancy_order.csv"))
summary <- list(source_interval_rows = nrow(p$data), source_pregnancy_rows = sum(p$data$pregnancy == 1, na.rm = TRUE),
  candidates_before_overlap = result$candidates_before_overlap,
  cleaned_episodes = nrow(e), pregnancy_eligible_women_before_TB_exclusions = nrow(result$pregnancy_eligible_observation),
  overlap_adjustments = result$overlap_adjustments, overlap_collapses = result$overlap_collapses,
  empty_after_observation_clipping = sum(result$audit$pregnancy_status == "empty_after_observation_clipping"),
  originally_left_truncated = sum(e$left_truncated), originally_right_truncated = sum(e$right_truncated),
  one_day_clipped_episodes = sum(e$observed_pregnancy_days == 1L),
  total_observed_pregnancy_days = sum(e$observed_pregnancy_days),
  first_two_episodes_flagged = sum(e$first_two_documented),
  note = "All cleaned pregnancy orders are retained. Prior TB exclusions, TB event ascertainment, and primary risk-set restrictions are later pieces.")
manifest <- finish_stage(ctx, inputs, outputs, rules, summary)
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Piece 3 complete: ", nrow(e), " cleaned episodes among ", nrow(result$pregnancy_eligible_observation), " women before prior-TB exclusions.")
