# Baseline covariates for the cohort models.
derive_baseline_covariates <- function(visits, cohort) {
  v <- data.table::copy(data.table::as.data.table(visits)[patient_id %in% cohort$patient_id,
    .(patient_id, source_row, apptdate, cd4date, tbProphStartDate, CLINIC, Weight, Height, bmi, cd4, WHOstage, viralload, Married, civilstatus)])
  v[, entry := cohort$entry[match(patient_id, cohort$patient_id)]]
  valid <- function(x, lo, hi) ifelse(is.finite(x) & x >= lo & x <= hi, x, NA_real_)
  v[, bmi_valid := valid(bmi, 10, 60)]
  v[is.na(bmi_valid), bmi_valid := valid(Weight / (Height / 100)^2, 10, 60)]
  v[, `:=`(cd4_valid = valid(cd4, 0, 5000), who_valid = valid(WHOstage, 1, 4),
    vl_valid = valid(viralload, 0, 100000000), married_derived = ifelse(Married %in% c(0, 1), Married, NA_real_))]
  v[!(WHOstage %in% 1:4), who_valid := NA_real_]
  v[is.na(married_derived), married_derived := c(0, 1, 1, 0, 0, 0)[match(civilstatus, 1:6)]]
  latest <- function(value_column, date_column, lookback = Inf) {
    z <- data.table::data.table(patient_id = v$patient_id, source_row = v$source_row,
      date = v[[date_column]], value = v[[value_column]], entry = v$entry)
    z <- z[!is.na(value) & !is.na(date) & date <= entry & as.numeric(entry - date) <= lookback]
    data.table::setorderv(z, c("patient_id", "date", "source_row"))
    z[, .SD[.N], by = patient_id]
  }
  specifications <- list(cd4 = c("cd4_valid", "cd4date"), who = c("who_valid", "apptdate"),
    bmi = c("bmi_valid", "apptdate"), vl = c("vl_valid", "apptdate"),
    clinic = c("CLINIC", "apptdate"), marital = c("married_derived", "apptdate"))
  out <- data.table::data.table(patient_id = cohort$patient_id)
  for (nm in names(specifications)) {
    spec <- specifications[[nm]]
    z <- latest(spec[1], spec[2], if (nm %in% c("clinic", "marital")) Inf else 365)
    idx <- match(out$patient_id, z$patient_id)
    data.table::set(out, j = paste0(nm, "_value"), value = z$value[idx])
    data.table::set(out, j = paste0(nm, "_date"), value = z$date[idx])
    data.table::set(out, j = paste0(nm, "_source_row"), value = z$source_row[idx])
  }
  common <- out[!is.na(clinic_value), .N, by = .(clinic = as.integer(clinic_value))][N >= 100, clinic]
  tpt_ids <- unique(v[!is.na(tbProphStartDate) & tbProphStartDate <= entry, patient_id])
  out[, `:=`(who_code = data.table::fcase(is.na(who_value), 2L, who_value >= 3, 1L, default = 0L),
    cd4_code = data.table::fcase(is.na(cd4_value), 4L, cd4_value < 200, 3L, cd4_value < 350, 2L, cd4_value < 500, 1L, default = 0L),
    bmi_code = data.table::fcase(is.na(bmi_value), 3L, bmi_value < 18.5, 1L, bmi_value >= 25, 2L, default = 0L),
    clinic_code = ifelse(as.integer(clinic_value) %in% common, as.integer(clinic_value), -1L),
    marital_code = data.table::fcase(is.na(marital_value), 2L, marital_value == 1, 1L, default = 0L),
    tpt_before_entry = as.integer(patient_id %in% tpt_ids))]
  data.table::setorderv(out, "patient_id")
  out
}
