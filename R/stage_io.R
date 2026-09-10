# Read and save derived datasets and check input file hashes.
start_stage <- function(project, stage) {
  for (p in c("data.table", "jsonlite", "digest"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing R package: ", p)
  data.table::setDTthreads(2L)
  Sys.setenv(TZ = "UTC")
  work <- file.path(project, "work", stage)
  report <- file.path(project, "outputs", stage)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  dir.create(report, recursive = TRUE, showWarnings = FALSE)
  list(project = project, stage = stage, work = work, report = report,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

checked_input <- function(project, stage, filename, expected) {
  path <- file.path(project, "work", stage, filename)
  if (is.null(expected) || !is.character(expected) || length(expected) != 1L)
    stop("Missing expected input checksum: ", stage, "/", filename)
  before <- digest::digest(path, file = TRUE, algo = "sha256")
  if (!identical(before, expected)) stop("Input checksum changed: ", stage, "/", filename, "; rerun its producing stage and validation.")
  list(data = readRDS(path), fingerprint = list(path = path, sha256 = before))
}

save_stage_rds <- function(data, directory, filename) {
  target <- file.path(directory, filename)
  temp <- paste0(target, ".partial")
  saveRDS(data, temp, compress = "gzip", version = 3)
  if (!file.rename(temp, target)) stop("Could not finalize derivative: ", filename)
  digest::digest(target, file = TRUE, algo = "sha256")
}

finish_stage <- function(ctx, input_fingerprints, output_hashes, rules, summary,
                         extra_scripts = character(), flags = list()) {
  for (x in input_fingerprints)
    if (!identical(digest::digest(x$path, file = TRUE, algo = "sha256"), x$sha256))
      stop("A standardized input changed during this stage.")
  args <- commandArgs(FALSE)
  script <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
  scripts <- c(normalizePath(script, mustWork = TRUE),
    file.path(ctx$project, "R", "eligibility_pregnancy.R"), file.path(ctx$project, "R", "stage_io.R"),
    file.path(ctx$project, extra_scripts))
  manifest <- list(stage = ctx$stage, started_utc = ctx$started_utc,
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    inputs = input_fingerprints, output_sha256 = output_hashes,
    code_sha256 = setNames(lapply(scripts, digest::digest, file = TRUE, algo = "sha256"), basename(scripts)),
    rules = rules, summary = summary, inputs_unchanged = TRUE,
    prior_TB_exclusions_applied = FALSE, TB_events_ascertained = FALSE,
    first_two_restriction_applied = FALSE, models_fitted = FALSE,
    R_version = R.version.string,
    package_versions = setNames(lapply(c("data.table", "jsonlite", "digest"),
      function(p) as.character(utils::packageVersion(p))), c("data.table", "jsonlite", "digest")))
  for (nm in names(flags)) manifest[[nm]] <- flags[[nm]]
  jsonlite::write_json(manifest, file.path(ctx$report, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null", digits = NA)
  writeLines(capture.output(sessionInfo()), file.path(ctx$report, "sessionInfo.txt"))
  invisible(manifest)
}
