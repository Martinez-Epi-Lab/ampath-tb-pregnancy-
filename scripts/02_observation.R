#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
source(file.path(project, "R", "stage_io.R"))
source(file.path(project, "R", "eligibility_pregnancy.R"))
ctx <- start_stage(project, "piece2")
prior <- jsonlite::read_json(file.path(project, "outputs", "piece1", "run_manifest.json"))
message("Piece 2: checking standardized inputs and deriving visit bounds.")
d <- checked_input(project, "piece1", "demographics.rds", prior$standardized_sha256$demographics)
v <- checked_input(project, "piece1", "visits.rds", prior$standardized_sha256$visits)
rules <- study_rules()
result <- observation_windows(d$data, v$data, rules)
inputs <- list(demographics = d$fingerprint, visits = v$fingerprint)
rm(v); invisible(gc())
outputs <- list()
for (nm in names(result)) outputs[[paste0(nm, ".rds")]] <- save_stage_rds(result[[nm]], ctx$work, paste0(nm, ".rds"))
flow <- result$audit[, .(people = .N, female_coded = sum(female_coded),
  male_coded = sum(!is.na(male) & male == 1), other_or_missing_sex = sum(is.na(male) | !male %in% c(0, 1))), by = observation_status]
setorderv(flow, "observation_status")
fwrite(flow, file.path(ctx$report, "eligibility_flow.csv"))
o <- result$observation
summary <- list(source_demographic_people = nrow(d$data), eligible_people_all_sexes = nrow(o),
  eligible_female_people = nrow(result$female_observation), eligible_male_people = sum(o$male == 1, na.rm = TRUE),
  one_day_observations = sum(o$observation_days == 1L), total_observation_days = sum(o$observation_days),
  earliest_entry = format(min(o$entry)), latest_exit = format(max(o$exit)),
  note = "The original observation algorithm precedes sex/pregnancy selection; its 90,926-person total is not a women-only count.")
manifest <- finish_stage(ctx, inputs, outputs, rules, summary)
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Piece 2 complete: ", nrow(o), " people with eligible observation, including ", nrow(result$female_observation), " women.")
