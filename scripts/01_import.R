#!/usr/bin/env Rscript
# Run from any working directory: Rscript --vanilla /path/to/scripts/01_import.R [config.R]
args_full <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_full[grepl("^--file=", args_full)][1])
this_file <- gsub("~+~", " ", this_file, fixed = TRUE)
project <- dirname(dirname(normalizePath(this_file, mustWork = TRUE)))
args <- commandArgs(trailingOnly = TRUE)
config <- if (length(args)) args[[1]] else file.path(project, "config", "paths.local.R")
if (!file.exists(config)) stop("Supply config/paths.local.R; see config/paths.example.R.")
needed <- c("haven", "data.table", "jsonlite", "digest")
for (pkg in needed) if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(2L)
Sys.setenv(TZ = "UTC")
source(file.path(project, "R", "import_helpers.R"))
cfg <- new.env(parent = baseenv())
sys.source(config, envir = cfg)
paths <- cfg$ampath_paths
if (!identical(sort(names(paths)), sort(names(import_schema())))) stop("Config must name the three core datasets.")
paths <- lapply(paths, normalizePath, mustWork = TRUE)
work <- file.path(project, "work", "piece1")
reports <- file.path(project, "outputs", "piece1")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
dir.create(reports, recursive = TRUE, showWarnings = FALSE)
if (any(vapply(paths, function(p) startsWith(p, paste0(normalizePath(work), "/")) ||
  startsWith(p, paste0(normalizePath(reports), "/")), logical(1)))) stop("Inputs cannot be stored in output directories.")
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
message("Calculating source-file SHA-256 hashes.")
before_hash <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
datasets <- list(); inventories <- list(); schemas <- list(); output_hashes <- list()
schema <- import_schema()
for (nm in names(schema)) {
  message("Importing ", nm, ".")
  # Read CSV fields as text before assigning their types.
  raw <- if (nm == "visits") data.table::fread(paths[[nm]], colClasses = "character",
    na.strings = NULL, strip.white = FALSE, check.names = FALSE, showProgress = FALSE) else
    data.table::as.data.table(haven::read_sas(paths[[nm]], encoding = "UTF-8", .name_repair = "minimal"))
  std <- standardize_source(raw, nm, schema[[nm]])
  inventories[[nm]] <- source_inventory(raw, std$data, nm)
  schemas[[nm]] <- std$schema
  datasets[[nm]] <- std$data[, .(patient_id)]
  target <- file.path(work, paste0(nm, ".rds"))
  temp <- paste0(target, ".partial")
  saveRDS(std$data, temp, compress = "gzip", version = 3)
  if (!file.rename(temp, target)) stop("Could not finalize standardized dataset.")
  output_hashes[[nm]] <- digest::digest(target, algo = "sha256", file = TRUE)
  message(nm, ": ", nrow(std$data), " rows retained.")
  rm(raw, std); invisible(gc())
}
inventory <- rbindlist(inventories)
column_schema <- rbindlist(schemas)
linkage <- linkage_inventory(datasets)
message("Checking source files remained byte-for-byte unchanged.")
after_hash <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
stopifnot(identical(before_hash, after_hash))
fwrite(inventory, file.path(reports, "dataset_inventory.csv"))
fwrite(column_schema, file.path(reports, "column_audit.csv"), na = "")
fwrite(linkage, file.path(reports, "linkage_audit.csv"))
manifest <- list(stage = "Step 1 import and standardization only", started_utc = started,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_files = lapply(names(paths), function(nm) list(dataset = nm, path = paths[[nm]],
    bytes = unname(file.info(paths[[nm]])$size), sha256_before = unname(before_hash[nm]),
    sha256_after = unname(after_hash[nm]))), standardized_sha256 = output_hashes,
  R_version = R.version.string, platform = R.version$platform,
  package_versions = setNames(lapply(needed, function(p) as.character(packageVersion(p))), needed),
  all_rows_retained = TRUE, source_files_unchanged = TRUE,
  no_eligibility_or_domain_filters_applied = TRUE,
  nonmissing_to_missing_cells = sum(column_schema$nonmissing_to_missing),
  notes = c("Participant-level standardized RDS files are local restricted intermediates.",
    "All original source columns retained; source_row maps to original file row order.",
    "No patient IDs or individual values are included in aggregate validation tables.",
    "Calendar-invalid dates become NA with original date field retained in source_*.",
    "Plausibility checks, cohort eligibility, pregnancy cleaning, and TB definitions are applied after import."))
jsonlite::write_json(manifest, file.path(reports, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(capture.output(sessionInfo()), file.path(reports, "sessionInfo.txt"))
writeLines(c("# Step 1 import summary", "", paste("Completed:", manifest$completed_utc), "",
  "All source rows were retained. No cohort eligibility, TB, or pregnancy-cleaning rules were applied.",
  "Source SHA-256 hashes matched before and after import.", "",
  paste0("- ", inventory$dataset, ": ", format(inventory$rows, big.mark = ",", trim = TRUE),
    " rows; ", inventory$unique_ids, " unique nonmissing IDs; ", inventory$duplicate_source_rows,
    " duplicate source rows retained."), "",
  paste("Nonmissing cells converted to missing:", manifest$nonmissing_to_missing_cells),
  "See column_audit.csv for exact per-column conversion and missingness counts and linkage_audit.csv for unmatched source records."),
  file.path(reports, "import_summary.md"))
message("Step 1 import complete. Aggregate reports: ", reports)
