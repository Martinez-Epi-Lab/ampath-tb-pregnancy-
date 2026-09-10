#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
this_file <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "stage2_io.R", "manuscript_rng.R", "pathway_biomarkers.R", "pathway_contacts.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage3_contacts"); check_core_versions(project); load_manuscript_rng(project)
c <- stage2_input(ctx, "piece4", "cohort.rds"); p <- stage2_input(ctx, "piece4", "episodes.rds")
audit <- stage2_input(ctx, "piece3", "audit.rds"); all_episodes <- stage2_input(ctx, "piece3", "episodes.rds")
day <- stage2_input(ctx, "stage3_contact_cells", "contact_dates.rds")
cells <- stage2_input(ctx, "stage3_contact_cells", "cells.rds")
message("Constructing within-episode contact periods.")
patient_order <- python_patient_order(audit$data, all_episodes$data)
phases <- build_pathway_episodes(c$data, p$data, patient_order)
phases <- phases[phase %in% c("preconception_days_180_to_1", "pregnancy", "postpartum_days_1_180")]
setnames(phases, c("phase_start", "phase_end"), c("start", "end"))
phases <- contact_episode_counts(phases, day$data, c$data)
pairs <- contact_pair_frames(phases)
rng <- new_manuscript_rng(20260829L); draws <- summaries <- list()
for (nm in names(pairs)) {
  message("Bootstrapping contact changes: ", nm, ".")
  draws[[nm]] <- bootstrap_contact_change(pairs[[nm]], rng)
  summaries[[nm]] <- summarize_contact_pairs(pairs[[nm]], draws[[nm]], nm)
}
message("Bootstrapping overall contact rates.")
draws$phase_rates <- bootstrap_contact_rates(cells$data, rng)
phase_summary <- contact_phase_summary(cells$data)
ci <- apply(draws$phase_rates, 2L, stats::quantile, c(0.025, 0.975), type = 7)
phase_summary[, `:=`(routine_visit_rate_ci_lower = ci[1, match(exposure, colnames(ci))],
  routine_visit_rate_ci_upper = ci[2, match(exposure, colnames(ci))])]
fwrite(phase_summary, file.path(ctx$report, "healthcare_contact_phase_rates.csv"))
fwrite(rbindlist(summaries), file.path(ctx$report, "healthcare_contact_paired_changes.csv"))
outputs <- list(episode_phase.rds = save_stage_rds(phases, ctx$work, "episode_phase.rds"),
  pairs.rds = save_stage_rds(pairs, ctx$work, "pairs.rds"), bootstrap_draws.rds = save_stage_rds(draws, ctx$work, "bootstrap_draws.rds"))
summary <- list(cohort_women = nrow(c$data), episodes = uniqueN(phases$episode_number),
  episodes_conceived_after_entry = uniqueN(phases[conceived_after_entry == TRUE, episode_number]),
  bootstrap_replicates = 5000L, seed = 20260829L)
finish_stage(ctx, list(cohort = c$fingerprint, episodes = p$fingerprint, pregnancy_audit = audit$fingerprint,
  all_episodes = all_episodes$fingerprint, contact_dates = day$fingerprint, cells = cells$fingerprint), outputs,
  list(bootstrap_unit = "woman", bootstrap_method = "PCG64 multinomial cluster counts; single continuous stream",
    replicates = 5000L, seed = 20260829L, quantile_type = 7L), summary,
  c("R/manuscript_rng.R", "src/manuscript_rng.c", "R/pathway_biomarkers.R", "R/pathway_contacts.R"),
  list(prior_TB_exclusions_applied = TRUE, TB_events_ascertained = TRUE, first_two_restriction_applied = TRUE))
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Contact rate and paired-change analyses complete.")
