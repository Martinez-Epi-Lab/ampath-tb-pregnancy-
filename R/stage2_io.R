# Input checks and result exports for secondary analyses.
stage2_input <- function(ctx, stage, filename) {
  manifest <- jsonlite::read_json(file.path(ctx$project, "outputs", stage, "run_manifest.json"))
  hash <- if (stage == "piece1") manifest$standardized_sha256[[sub("[.]rds$", "", filename)]] else manifest$output_sha256[[filename]]
  checked_input(ctx$project, stage, filename, hash)
}

check_core_versions <- function(project) {
  m <- lapply(c("piece1", "piece4", "piece5", "piece6"), function(s) jsonlite::read_json(file.path(project, "outputs", s, "run_manifest.json")))
  stopifnot(identical(m[[2]]$inputs$visits$sha256, m[[1]]$standardized_sha256$visits),
    identical(m[[3]]$inputs$cohort$sha256, m[[2]]$output_sha256$cohort.rds),
    identical(m[[3]]$inputs$episodes$sha256, m[[2]]$output_sha256$episodes.rds),
    identical(m[[4]]$inputs$cohort$sha256, m[[2]]$output_sha256$cohort.rds),
    identical(m[[4]]$inputs$visits$sha256, m[[1]]$standardized_sha256$visits))
}

stage2_model_exports <- function(ctx, models, descriptions) {
  coefficients <- data.table::rbindlist(lapply(names(models), function(id) {
    x <- data.table::copy(models[[id]]$estimates); x[, model_id := id]; x
  }), fill = TRUE)
  diagnostics <- data.table::rbindlist(lapply(names(models), function(id) c(list(model_id = id), models[[id]]$diagnostics)), fill = TRUE)
  data.table::fwrite(coefficients, file.path(ctx$report, "model_coefficients.csv"))
  data.table::fwrite(diagnostics, file.path(ctx$report, "model_diagnostics.csv"))
  if (length(descriptions)) {
    summaries <- data.table::rbindlist(descriptions, fill = TRUE)
    summaries <- merge(summaries, coefficients, by = c("model_id", "feature"), all.x = TRUE, sort = FALSE)
    summaries[risk_code == 0, irr := 1]
    summaries[, estimate_display := ifelse(risk_code == 0, "1.00 (Reference)", sprintf("%.2f (%.2f–%.2f)", irr, ci_lower, ci_upper))]
    data.table::fwrite(summaries, file.path(ctx$report, "exposure_estimates.csv"))
  }
  for (id in names(models)) {
    v <- models[[id]]$covariance
    data.table::fwrite(data.table::data.table(feature = rownames(v), data.table::as.data.table(v)), file.path(ctx$report, paste0(id, "_covariance.csv")))
  }
}
