# Analysis modules in dependency order.
pipeline_steps <- function() {
  data.frame(
    number = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
    name = c("import", "observation", "pregnancies", "tb", "person_time", "primary_models", "pregnancy_outcomes", "secondary_models", "postpartum_spline", "biomarkers", "contact_cells", "contact_rates", "mediation", "characteristics", "tables", "figures", "documents", "result_dictionary"),
    script = c("scripts/01_import.R", "scripts/02_observation.R", "scripts/03_pregnancies.R", "scripts/04_tb.R", "scripts/05_person_time.R", "scripts/06_primary_models.R", "scripts/07_pregnancy_outcomes.R", "scripts/08_secondary_models.R", "scripts/09_postpartum_spline.R", "scripts/10_biomarkers.R", "scripts/11_contact_cells.R", "scripts/12_contact_rates.R", "scripts/13_mediation.R", "scripts/14_characteristics.R", "scripts/15_tables.R", "scripts/16_figures.R", "scripts/17_documents.R", "scripts/18_result_dictionary.R"),
    description = c("Import source extracts", "Define observation", "Clean pregnancy episodes", "Ascertain TB", "Construct person-time", "Fit primary models", "Link pregnancy outcomes", "Fit secondary models", "Fit postpartum splines", "Analyze biomarker changes", "Construct contact cells", "Analyze contact rates", "Fit mediation models", "Summarize characteristics", "Build tables", "Draw figures", "Write editable tables", "Export result dictionaries"),
    stringsAsFactors = FALSE
  )
}
