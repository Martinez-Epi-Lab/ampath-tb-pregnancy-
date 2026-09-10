#!/usr/bin/env Rscript
# Export field definitions, source checksums, and numerical results.
args <- commandArgs(FALSE)
self <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(self, mustWork = TRUE)))
suppressPackageStartupMessages(library(data.table))
report <- file.path(project, "outputs/stage4_publication")
dir.create(report, recursive = TRUE, showWarnings = FALSE)
audit <- fread(file.path(project, "outputs/piece1/column_audit.csv"))
dictionary <- audit[, .(dataset, column, required, intended_type, standardized_class, source_label, source_format)]
outcome <- data.table(dataset = 'pregnancy_outcomes', column = c('patient_id','apptdate','pregnancyoutcome'),
  required = TRUE, intended_type = c('identifier','date','numeric'), standardized_class = c('character','Date','integer'),
  source_label = c('Patient/date row keys must match the visits extract','Appointment date; explicit day-month-year parsing','Recorded outcome codes 1–8; see R/pregnancy_outcomes.R for adjudication'), source_format = '')
fwrite(rbindlist(list(dictionary,outcome)), file.path(report,'source_data_dictionary.csv'))
source(file.path(project,'R/data_dictionary.R'))
derived <- derived_variable_dictionary()
fwrite(derived,file.path(report,'derived_data_dictionary.csv'))
imports <- jsonlite::read_json(file.path(project,'outputs/piece1/run_manifest.json'))
source_hashes <- rbindlist(lapply(imports$source_files,function(x)data.table(dataset=x$dataset,filename=basename(x$path),sha256=x$sha256_before,unchanged=identical(x$sha256_before,x$sha256_after))))
o <- jsonlite::read_json(file.path(project,'outputs/stage2_outcomes/run_manifest.json'))
source_hashes <- rbind(source_hashes,data.table(dataset='pregnancy_outcomes',filename=basename(o$inputs$outcome_source$path),sha256=o$inputs$outcome_source$sha256,unchanged=o$summary$source_unchanged))
stopifnot(nrow(source_hashes)==4L,all(source_hashes$unchanged))
fwrite(source_hashes,file.path(report,'source_checksums.csv'))
software <- data.table(component=c('R','platform','haven','data.table','jsonlite','digest','survival'),
  version=c(as.character(getRversion()),R.version$platform,vapply(c('haven','data.table','jsonlite','digest','survival'),function(p)as.character(utils::packageVersion(p)),character(1))))
fwrite(software,file.path(report,'software_versions.csv'))
values <- fread(file.path(project,'outputs/stage4_tables/numerical_results_register.csv'))
extras <- c(table1='stage4_table1/table1_statistics.csv',flow='stage4_table1/flow_counts.csv',
  spline_curve='stage2_spline/postpartum_spline_curve.csv',spline_tests='stage2_spline/postpartum_spline_tests.csv',
  postpartum_tests='stage2_secondary/postpartum_interval_tests.csv',outcomes='stage2_outcomes/outcome_completeness.csv',
  baseline_missingness='stage4_table1/table1_missingness.csv',incidence='stage4_tables/incidence_intervals.csv')
for (key in names(extras)) {
  d <- fread(file.path(project,'outputs',extras[[key]]))
  for (col in names(d)[vapply(d,is.numeric,logical(1))]) {
    get <- function(nm) if (nm %in% names(d)) as.character(d[[nm]]) else rep('',nrow(d))
    values <- rbind(values,data.table(source=key,source_row=seq_len(nrow(d)),measure=col,value=d[[col]],
      model=get('model_id'),feature=if('characteristic'%in%names(d))paste(get('population'),get('characteristic'),get('level'),sep=' | ')else if('key'%in%names(d))get('key')else get('feature'),comparison=get('comparison'),phase=get('phase')))
  }
}
fwrite(values,file.path(report,'numerical_results_register.csv'))
source_map <- rbind(data.table(source=c('core','secondary','influence','coverage','biomarker_pairs','contact_rates','contact_pairs','mediation'),
  path=c('piece6/exposure_estimates.csv','stage2_secondary/exposure_estimates.csv','stage2_secondary/covariate_influence.csv','stage3_biomarkers/biomarker_phase_coverage.csv','stage3_biomarkers/biomarker_paired_changes.csv','stage3_contacts/healthcare_contact_phase_rates.csv','stage3_contacts/healthcare_contact_paired_changes.csv','stage3_mediation/mediation_effects.csv')),
  data.table(source=names(extras),path=unname(extras)))
fwrite(source_map,file.path(report,'result_source_map.csv'))
jsonlite::write_json(list(original_fields=nrow(dictionary),additional_outcome_fields=nrow(outcome),derived_definitions=nrow(derived),numerical_values=nrow(values),source_files=nrow(source_hashes)),file.path(report,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
message('Exported data dictionaries, numerical results, source checksums, and software versions.')
