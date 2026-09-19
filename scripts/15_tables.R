#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
self <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(self, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
for (f in c("stage_io.R", "manuscript_tables.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage4_tables")
inputs <- list(); dir.create(file.path(ctx$report, "full_precision"), showWarnings = FALSE)
read_csv <- function(stage, file, key) {
  path <- file.path(project, "outputs", stage, file)
  inputs[[key]] <<- list(path = path, sha256 = digest::digest(path, file = TRUE, algo = "sha256"))
  x <- fread(path); fwrite(x, file.path(ctx$report, "full_precision", paste0(key, ".csv"))); x
}
x <- list(core = read_csv("piece6", "exposure_estimates.csv", "core"),
  secondary = read_csv("stage2_secondary", "exposure_estimates.csv", "secondary"),
  influence = read_csv("stage2_secondary", "covariate_influence.csv", "influence"),
  table1 = read_csv("stage4_table1", "table1.csv", "table1"),
  coverage = read_csv("stage3_biomarkers", "biomarker_phase_coverage.csv", "biomarker_coverage"),
  biomarker_pairs = read_csv("stage3_biomarkers", "biomarker_paired_changes.csv", "biomarker_pairs"),
  contact_rates = read_csv("stage3_contacts", "healthcare_contact_phase_rates.csv", "contact_rates"),
  contact_pairs = read_csv("stage3_contacts", "healthcare_contact_paired_changes.csv", "contact_pairs"),
  mediation = read_csv("stage3_mediation", "mediation_effects.csv", "mediation"))
outcome <- file.path(project, "outputs/stage2_outcomes/summary.json")
inputs$outcome_summary <- list(path = outcome, sha256 = digest::digest(outcome, file = TRUE, algo = "sha256"))
x$outcome_summary <- jsonlite::read_json(outcome)
tables <- build_manuscript_tables(x)
outputs <- list(tables.rds = save_stage_rds(tables, ctx$work, "tables.rds"))
registry <- cells <- list()
for (id in names(tables)) {
  tab <- tables[[id]]; fwrite(tab$data, file.path(ctx$report, paste0(id, ".csv")))
  registry[[id]] <- data.table(table_id = id, title = tab$title, rows = nrow(tab$data), columns = ncol(tab$data), landscape = tab$landscape, notes = tab$notes)
  cells[[id]] <- rbindlist(lapply(seq_len(ncol(tab$data)), function(j) data.table(table_id = id, row = seq_len(nrow(tab$data)), column = j,
    column_label = names(tab$data)[j], value = tab$data[[j]])))
}
fwrite(rbindlist(registry), file.path(ctx$report, "table_register.csv"))
fwrite(rbindlist(cells), file.path(ctx$report, "table_cells.csv"))
# Collect numeric results with their source rows and model labels.
values <- list()
for (key in setdiff(names(x), c("table1", "outcome_summary"))) {
  data <- x[[key]]
  for (col in names(data)[vapply(data, is.numeric, logical(1))]) {
    index <- length(values) + 1L
    values[[index]] <- data.table(source = key, source_row = seq_len(nrow(data)), measure = col, value = data[[col]],
      model = if ("model_id" %in% names(data)) data$model_id else "",
      feature = if ("feature" %in% names(data)) data$feature else if ("exposure" %in% names(data)) data$exposure else "",
      comparison = if ("comparison" %in% names(data)) data$comparison else "", phase = if ("phase" %in% names(data)) data$phase else "")
  }
}
fwrite(rbindlist(values), file.path(ctx$report, "numerical_results_register.csv"))
incidence <- rbindlist(lapply(c("other_observed_time", "pregnancy", "postpartum_days_1_180"), function(f) {
  z <- x$core[model_id == "cohort_primary" & feature == f]; cbind(data.table(feature = f), poisson_incidence(z$events, z$person_days))
}))
incidence[, method := "Wilson–Hilferty approximation to the Poisson interval, matching the original calculation"]
fwrite(incidence, file.path(ctx$report, "incidence_intervals.csv"))
summary <- list(main_tables = 3L, supplementary_tables = 11L, physical_tables = length(tables), formatted_cells = nrow(rbindlist(cells)),
  numerical_values = nrow(rbindlist(values)))
finish_stage(ctx, inputs, outputs, list(rounding = "Decimal half up", table3_strict_tb_row = FALSE), summary, "R/manuscript_tables.R")
jsonlite::write_json(summary, file.path(ctx$report, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Prepared all 3 main and 11 supplementary tables.")
