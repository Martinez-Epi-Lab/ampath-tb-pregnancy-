# Derived variables

| dataset | variable | type | units_or_categories | definition | code |
|---|---|---|---|---|---|
| observation | patient_id | character | key | Canonical source identifier; values remain restricted | R/import_helpers.R |
| observation | entry | Date | date | First eligible observation day after applying run-in and age/calendar entry rules | R/eligibility_pregnancy.R |
| cohort | administrative_exit | Date | date | Last eligible observation day before TB censoring | R/tb_events.R |
| cohort | event_date | Date | date | Earliest qualifying primary TB date within eligible follow-up | R/tb_events.R |
| cohort | followup_end | Date | date | Earlier of administrative exit and first qualifying TB event | R/tb_events.R |
| pregnancies | original_start | Date | date | Pregnancy start supplied in the interval extract | R/eligibility_pregnancy.R |
| pregnancies | original_end | Date | date | Pregnancy end supplied in the interval extract | R/eligibility_pregnancy.R |
| pregnancies | analysis_start | Date | date | Later of cleaned episode start and observation entry | R/eligibility_pregnancy.R |
| pregnancies | analysis_end | Date | date | Earlier of overlap-truncated episode end and observation exit | R/eligibility_pregnancy.R |
| pregnancies | documented_pregnancy_order | integer | order | Order assigned after pregnancy eligibility and overlap processing | R/eligibility_pregnancy.R |
| person_time | days | integer | days | Inclusive observed interval duration | R/person_time.R |
| person_time | events | integer | count | Qualifying TB events allocated to the retained interval | R/person_time.R |
| person_time | risk_code | integer | model-specific exposure category | Exposure category defined by the selected pooled or windowed analysis | R/person_time.R |
| person_time | age_code | integer | age group | Attained age category used by model_matrix | R/models.R |
| person_time | calendar_code | integer | calendar group | Calendar period category used by model_matrix | R/models.R |
| baseline | who_value | numeric | recorded WHO stage | Most recent valid baseline WHO stage in the defined lookback window | R/baseline_covariates.R |
| baseline | cd4_value | numeric | cells per microliter | Most recent valid baseline CD4 count | R/baseline_covariates.R |
| baseline | bmi_value | numeric | kg per square meter | Most recent valid baseline BMI | R/baseline_covariates.R |
| baseline | marital_value | numeric | derived binary married status | Most recent valid Married indicator with civilstatus recoding as fallback | R/baseline_covariates.R |
| baseline | tpt_before_entry | integer | 0 or 1 | Tuberculosis preventive treatment recorded on or before study entry | R/baseline_covariates.R |
| outcomes | recorded_live_birth | logical | FALSE or TRUE | Episode has a recorded live birth under the primary outcome adjudication | R/pregnancy_outcomes.R |
| outcomes | primary90_vital_outcome | character | adjudicated category | Vital outcome using linked records from 30 days before through 90 days after original pregnancy end | R/pregnancy_outcomes.R |
| pathways | phase | character | episode phase | Observed preconception or pregnancy or postpartum window clipped to observation and TB censoring | R/pathway_biomarkers.R |
| contacts | routine_visit_dates | integer | count | Distinct contact dates without TB encounter or investigation entries and excluding the TB event day | R/pathway_contacts.R |
| contacts | all_visit_dates | integer | count | Distinct contact dates excluding the TB event day | R/pathway_contacts.R |
| mediation | mediator | numeric | log transformed contact rate | Log of 1 plus distinct qualifying contact dates per 100 observed days | R/pathway_mediation.R |
| results | irr | numeric | rate ratio | Exponentiated model coefficient or specified contrast | R/models.R |
| results | ci_lower | numeric | rate ratio | Lower 95 percent confidence limit | R/models.R |
| results | ci_upper | numeric | rate ratio | Upper 95 percent confidence limit | R/models.R |
| Table1 | pregnancy_count | integer | 0 or 1 or 2 | First two documented pregnancies contributing observed pregnancy time in the specified design | R/table1.R |
| Table1 | retained_analysis_years | numeric | years | Design-specific retained analysis days divided by 365.25 | R/table1.R |
