#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
source(file.path(project, "R", "stage_io.R"))
source(file.path(project, "R", "eligibility_pregnancy.R"))
source(file.path(project, "R", "person_time.R"))
ctx <- start_stage(project, "piece5")
p4 <- jsonlite::read_json(file.path(project, "outputs/piece4/run_manifest.json"))
c <- checked_input(project, "piece4", "cohort.rds", p4$output_sha256$cohort.rds)
p <- checked_input(project, "piece4", "episodes.rds", p4$output_sha256$episodes.rds)
message("Piece 5: constructing disjoint person-time intervals and model views.")
result <- construct_person_time(c$data, p$data)
outputs <- list(); descriptions <- list()
for (nm in names(result)) outputs[[paste0(nm, ".rds")]] <- save_stage_rds(result[[nm]], ctx$work, paste0(nm, ".rds"))
for (design in c("cohort", "sccs")) for (view in c("pooled_primary", "pooled_five", "postpartum_30day")) {
  rows <- analysis_rows(result[[paste0(design, "_intervals")]], view, sccs = design == "sccs")
  nm <- paste0(design, "_", view, ".rds")
  outputs[[nm]] <- save_stage_rds(rows, ctx$work, nm)
  descriptions[[nm]] <- describe_person_time(rows, design, view)
}
fwrite(rbindlist(descriptions, fill = TRUE), file.path(ctx$report, "person_time_summary.csv"))
fwrite(result$event_audit[, .(events = .N), by = .(risk_code, order_code, included_in_primary_analysis)], file.path(ctx$report, "event_flow.csv"))
summary <- list(cohort_women = nrow(c$data), total_first_TB_events = nrow(result$event_audit),
  analyzed_events = sum(result$event_audit$included_in_primary_analysis),
  excluded_later_pregnancy_events = sum(!result$event_audit$included_in_primary_analysis),
  cohort_person_days = result$cohort_intervals[included == TRUE, sum(days)],
  cohort_excluded_person_days = result$cohort_intervals[included == FALSE, sum(days)],
  sccs_person_days = result$sccs_intervals[included == TRUE & primary_case == TRUE, sum(days)],
  SCCS_cases = uniqueN(result$sccs_intervals[primary_case == TRUE, patient_id]))
finish_stage(ctx, list(cohort = c$fingerprint, episodes = p$fingerprint), outputs, person_time_rules(), summary,
  "R/person_time.R", list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, first_two_restriction_applied = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
message("Piece 5 complete: ", summary$analyzed_events, " analyzed events; ", summary$cohort_person_days, " cohort person-days.")
