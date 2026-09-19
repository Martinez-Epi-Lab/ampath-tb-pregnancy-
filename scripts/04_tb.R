#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
source(file.path(project, "R", "stage_io.R"))
source(file.path(project, "R", "eligibility_pregnancy.R"))
source(file.path(project, "R", "tb_events.R"))
ctx <- start_stage(project, "piece4")
p1 <- jsonlite::read_json(file.path(project, "outputs/piece1/run_manifest.json"))
p3 <- jsonlite::read_json(file.path(project, "outputs/piece3/run_manifest.json"))
d <- checked_input(project, "piece1", "demographics.rds", p1$standardized_sha256$demographics)
v <- checked_input(project, "piece1", "visits.rds", p1$standardized_sha256$visits)
o <- checked_input(project, "piece3", "pregnancy_eligible_observation.rds", p3$output_sha256$pregnancy_eligible_observation.rds)
p <- checked_input(project, "piece3", "episodes.rds", p3$output_sha256$episodes.rds)
p2 <- jsonlite::read_json(file.path(project, "outputs/piece2/run_manifest.json"))
stopifnot(identical(p3$inputs$demographics$sha256, d$fingerprint$sha256),
  identical(p2$inputs$visits$sha256, v$fingerprint$sha256),
  identical(p3$inputs$observation$sha256, p2$output_sha256$observation.rds))
message("Step 4: extracting TB evidence and corroborating pairs.")
evidence <- extract_tb_evidence(v$data)
tb <- ascertain_tb(evidence)
result <- construct_tb_cohort(o$data, d$data, p$data, tb$dates)
result$evidence <- evidence; result$tb_dates <- tb$dates; result$tb_candidates <- tb$candidates
outputs <- list()
for (nm in names(result)) outputs[[paste0(nm, ".rds")]] <- save_stage_rds(result[[nm]], ctx$work, paste0(nm, ".rds"))
fwrite(result$audit[, .(women = .N), by = cohort_status], file.path(ctx$report, "cohort_flow.csv"))
fwrite(evidence[, .(evidence_records = .N, people = uniqueN(patient_id)), by = source], file.path(ctx$report, "evidence_summary.csv"))
summary <- list(pregnancy_eligible_women = nrow(o$data), final_cohort_women = nrow(result$cohort),
  incident_first_events_all_orders = sum(!is.na(result$cohort$event_date)),
  cleaned_pregnancy_episodes_final_cohort = nrow(result$episodes),
  first_two_pregnancies_final_cohort = sum(result$episodes$first_two_documented),
  evidence_records = nrow(evidence), qualifying_candidate_records = nrow(tb$candidates))
finish_stage(ctx, list(demographics = d$fingerprint, visits = v$fingerprint, observation = o$fingerprint, episodes = p$fingerprint),
  outputs, tb_rules(), summary, "R/tb_events.R",
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
message("Step 4 complete: ", summary$final_cohort_women, " women; ", summary$incident_first_events_all_orders, " first TB events before the first-two restriction.")
