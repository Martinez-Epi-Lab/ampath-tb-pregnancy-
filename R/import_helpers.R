# Import source data before applying eligibility criteria.

import_schema <- function() {
  list(
    demographics = list(
      dates = c("dob", "dod", "arvstart", "arvstartTHERAP"),
      character = "pttype",
      numeric = c("birthdate_estimated", "Death", "presumdead", "lost2fup",
        "TravelTime", "AttendedSchool", "YearsOfSchool", "TimesPregnant",
        "CurrentlyInSchool", "male", "everpreg", "deathdate_estimated", "txnaive",
        "ageatarvstartTHERAP", "inhnumadmeas", "INHperfectad", "ageatenroll",
        "dropout", "dropoutreason", "tbtxprior", "ptidno", "ageatarvstart")),
    pregnancies = list(dates = c("intstart", "intend"), character = character(),
      numeric = c("estconMEDIAN", "estconTERM", "endEQdelp", "pregnancy", "intdur")),
    visits = list(
      dates = c("apptdate", "cd4date", "tbtxstartdate", "TBtxCompDate",
        "tbProphStartDate", "tbtxoutcomeDate"),
      character = c("cxr", "combined_enc_types", "tbextrapulmonaryPC",
        "tbpulmonaryPC", "ZNstain", "AFBculture", "regimen"),
      numeric = c("CLINIC", "Weight", "Height", "AlcoholHowOften",
        "AlcoholNumberOfDrinks", "AlcoholSixOrMoreDrinks", "bmi", "TBprophstop",
        "pmtctTX", "cd4per", "cd4", "WHOstage", "onARV", "tbtx",
        "tbextrapulmdiag", "tbpulmdiag", "hepatitis", "TBTreatmentCompleted",
        "TBtxCompletedYrs", "Married", "age", "encounter_type", "diabetes",
        "cxrNormal", "cxrInfiltrate", "cxrCavitary", "HbA1c", "missedTBprophLastMon",
        "inhpillstakenlast7days", "drinkalcohol", "cigsmoke", "cignumber", "cigdur",
        "cigstop", "weightgain", "TBprophadherence", "viralload", "VLundetectable",
        "lymphnodeTB", "lymphnodeTBDiag", "momTBtx", "careStatus", "fupstatus",
        "currentEDUC", "artcategory", "tbtxoutcome", "typetbpt", "tbptref",
        "visittype", "civilstatus", "transfer", "transferin", "catchmentarea",
        "incomesource", "tbpulmonary", "tbextrapulmonary", "activeTB", "tbproph",
        "inh", "_3HP"))
  )
}

# Recognize CSV missing-value tokens and record SAS numeric missing codes separately.
# Keep clinical codes such as 0, 9, and 99.
pandas_na_tokens <- function() c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND",
  "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA",
  "NULL", "NaN", "None", "n/a", "nan", "null")

blank_to_na <- function(x) {
  if (!is.character(x)) return(x)
  x[x %in% pandas_na_tokens()] <- NA_character_
  x[!is.na(x) & grepl("^[[:space:]]*$", x)] <- NA_character_
  x
}

canonical_id <- function(x) {
  if (is.numeric(x)) {
    if (any(!is.na(x) & (!is.finite(x) | x < 0 | x != floor(x) | x > 2^53 - 1)))
      stop("Invalid numeric patient_id; must be a nonnegative exactly representable integer.")
    y <- rep(NA_character_, length(x))
    y[!is.na(x)] <- sprintf("%.0f", x[!is.na(x)])
    return(y)
  }
  y <- blank_to_na(trimws(as.character(x)))
  # Accept IDs ending in .0 without converting them to floating-point numbers.
  good <- is.na(y) | grepl("^[0-9]+(?:\\.0+)?$", y, perl = TRUE)
  if (any(!good)) stop("Invalid character patient_id; nonintegral or unexpected syntax.")
  y <- sub("\\.0+$", "", y)
  y <- sub("^0+(?=[0-9])", "", y, perl = TRUE)
  y
}

parse_ampath_date <- function(x) {
  # Parse English month names regardless of the computer's language settings.
  raw <- blank_to_na(as.character(x))
  unique_raw <- unique(raw[!is.na(raw)])
  y <- toupper(unique_raw)
  valid <- grepl("^[0-9]{1,2}[A-Z]{3}[0-9]{4}$", y)
  iso <- rep(NA_character_, length(y))
  if (any(valid)) {
    z <- y[valid]
    day <- as.integer(sub("([0-9]{1,2}).*", "\\1", z))
    mon <- match(sub("^[0-9]{1,2}([A-Z]{3}).*", "\\1", z), toupper(month.abb))
    yr <- as.integer(sub(".*([0-9]{4})$", "\\1", z))
    iso[valid] <- sprintf("%04d-%02d-%02d", yr, mon, day)
  }
  parsed <- as.Date(iso, format = "%Y-%m-%d")
  # Reject invalid dates, including dates automatically shifted into another month.
  bad <- !is.na(parsed) & format(parsed, "%Y-%m-%d") != iso
  parsed[bad] <- as.Date(NA)
  parsed[match(raw, unique_raw)]
}

normalize_sas_date <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  if (inherits(x, "POSIXt")) {
    out <- as.Date(x, tz = "UTC")
    if (any(!is.na(x) & as.numeric(x) != as.numeric(as.POSIXct(out, tz = "UTC"))))
      stop("SAS datetime contains nonmidnight times; needs an explicit conversion rule.")
    return(out)
  }
  # Use SAS date metadata to distinguish dates from numeric values.
  stop("Expected a SAS Date or POSIXt column; unrecognized date storage.")
}

parse_numeric_field <- function(x) {
  if (is.numeric(x)) return(haven::zap_missing(as.numeric(x)))
  y <- blank_to_na(as.character(x))
  # Treat SAS special missing codes as missing values.
  y[!is.na(y) & grepl("^\\.[A-Za-z_]?$", y)] <- NA_character_
  z <- suppressWarnings(as.numeric(y))
  z[!is.finite(z)] <- NA_real_
  z
}

date_range <- function(x) {
  if (all(is.na(x))) return(c(NA_character_, NA_character_))
  format(range(x, na.rm = TRUE), "%Y-%m-%d")
}

standardize_source <- function(raw, dataset, schema) {
  required <- c("patient_id", schema$dates, schema$character, schema$numeric)
  if (anyDuplicated(names(raw))) stop(dataset, ": duplicated source column names.")
  absent <- setdiff(required, names(raw))
  if (length(absent)) stop(dataset, ": missing required columns: ", paste(absent, collapse = ", "))
  if (any(grepl("^source_", names(raw)))) stop("Source names collide with provenance columns.")
  source_names <- names(raw)
  out <- data.table::copy(raw)
  data.table::set(out, j = "source_row", value = seq_len(nrow(raw)))
  data.table::set(out, j = "source_patient_id", value = raw[["patient_id"]])
  data.table::set(out, j = "patient_id", value = canonical_id(raw[["patient_id"]]))
  rows <- list()
  for (nm in source_names) {
    before <- raw[[nm]]
    type <- if (nm == "patient_id") "identifier" else if (nm %in% schema$dates) "date" else
      if (nm %in% schema$numeric) "numeric" else "character"
    source_missing <- if (is.character(before)) is.na(blank_to_na(before)) else is.na(before)
    tagged <- if (is.double(before)) sum(haven::is_tagged_na(before)) else 0L
    if (type == "date") {
      data.table::set(out, j = paste0("source_", nm), value = before)
      value <- if (dataset == "visits") parse_ampath_date(before) else normalize_sas_date(before)
      data.table::set(out, j = nm, value = value)
    } else if (type == "numeric") {
      value <- parse_numeric_field(before)
      data.table::set(out, j = nm, value = value)
    } else if (type == "character") {
      data.table::set(out, j = nm, value = blank_to_na(as.character(before)))
    }
    after <- out[[nm]]
    failures <- sum(!source_missing & is.na(after))
    # Keep the source value when conversion produces NA.
    if (failures && type != "date") data.table::set(out, j = paste0("source_", nm), value = before)
    rng <- if (type == "date") date_range(after) else c(NA_character_, NA_character_)
    rows[[length(rows) + 1L]] <- data.table::data.table(dataset = dataset, column = nm,
      required = nm %in% required, intended_type = type,
      source_class = paste(class(before), collapse = "/"),
      standardized_class = paste(class(after), collapse = "/"),
      source_label = if (is.null(attr(before, "label"))) NA_character_ else attr(before, "label"),
      source_format = if (is.null(attr(before, "format.sas"))) NA_character_ else attr(before, "format.sas"),
      rows = length(before), source_missing = sum(source_missing),
      standardized_missing = sum(is.na(after)), nonmissing_to_missing = failures,
      sas_tagged_missing = tagged, date_min = rng[1], date_max = rng[2])
  }
  stopifnot(nrow(raw) == nrow(out), identical(out$source_row, seq_len(nrow(raw))))
  list(data = out, schema = data.table::rbindlist(rows))
}

source_inventory <- function(raw, standardized, dataset) {
  id <- standardized$patient_id
  # Check source records before adding row numbers or changing field types.
  data.table::data.table(dataset = dataset, rows = nrow(raw), source_columns = ncol(raw),
    standardized_columns = ncol(standardized), unique_ids = data.table::uniqueN(id, na.rm = TRUE),
    missing_ids = sum(is.na(id)), duplicate_source_rows = sum(duplicated(raw)),
    repeated_id_rows = sum(duplicated(id[!is.na(id)])))
}

linkage_inventory <- function(datasets) {
  ans <- list()
  for (from in names(datasets)) for (to in setdiff(names(datasets), from)) {
    x <- datasets[[from]]$patient_id
    ref <- unique(datasets[[to]]$patient_id)
    ref <- ref[!is.na(ref)]
    matched <- !is.na(x) & x %in% ref
    ans[[length(ans) + 1L]] <- data.table::data.table(from = from, to = to,
      source_rows = length(x), matched_rows = sum(matched),
      unmatched_nonmissing_rows = sum(!is.na(x) & !matched), missing_id_rows = sum(is.na(x)),
      matched_unique_ids = data.table::uniqueN(x[matched]),
      unmatched_unique_ids = data.table::uniqueN(x[!is.na(x) & !matched]))
  }
  data.table::rbindlist(ans)
}
