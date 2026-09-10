# Analysis methods and code map

## Observation, pregnancy, and TB

The import preserves source rows and fields, standardizes identifiers and dates, and records missingness and conversion checks. Source clinical codes are not globally recoded as missing. Pregnancy-interval duration in the supplied SAS extract is measured in weeks. The field contract is defined in `R/import_helpers.R`.

`R/eligibility_pregnancy.R` defines observation within January 1, 2010–April 15, 2025, with a 90-day run-in and age limits of 15 and 50 years. It excludes invalid pregnancy intervals, truncates an earlier overlapping episode before the next pregnancy, removes episodes that collapse, assigns documented pregnancy order, and clips retained intervals to observation.

`R/tb_events.R` combines treatment starts, active-TB and treatment reports, diagnoses, bacteriological evidence, and corroborating records. Previous treatment and broad evidence before entry determine exclusions. Primary and stricter TB definitions remain separate.

`R/person_time.R` uses inclusive days. Primary analyses retain the first two documented pregnancies, postpartum days 1–180, and other eligible observation. Postpartum ends before a subsequent pregnancy. Cohort follow-up stops at the first TB event; SCCS retains full observation for cases whose event occurs in retained time. An event registry can therefore contain events excluded by a particular analysis restriction.

## Baseline characteristics and models

`R/baseline_covariates.R` selects the latest eligible baseline measurements. WHO stage, CD4, BMI, and viral load use the specified one-year lookback; clinic and marital status use available pre-entry history. Missing categories remain explicit. Preventive treatment recorded on entry is included.

`R/models.R` fits Poisson cohort models and conditional Poisson SCCS models with log-person-days offsets and model-specific age, calendar, and clinical covariates. Cohort covariance uses a woman-clustered sandwich estimator; conditional SCCS uses inverse-information covariance.

`R/secondary_analysis.R` implements expanded adjustment, covariate influence, all documented pregnancies, recorded live birth, pregnancy order, ART timing/duration, fully observed pregnancy, postpartum lags and 60-day windows, preconception exclusions, calendar restrictions, and stricter TB. The stricter-TB result remains available in the analysis output but is omitted from the manuscript Table 3 layout.

`R/pregnancy_outcomes.R` requires exact patient/date row alignment between the visits and outcome extracts. Records link to original pregnancy ends using deterministic nearest-episode rules. Primary adjudication uses days −30 through +90; extended linkage reaches +180. Missing outcome documentation is not interpreted as live birth.

`R/postpartum_spline.R` defines the restricted cubic spline basis, observed person-time, daily predictions, confidence intervals, and Wald tests.

## Pathway analyses

`R/pathway_biomarkers.R` constructs episode-specific preconception, pregnancy, and postpartum windows clipped to observation and first TB censoring. It selects the latest valid phase measurement and calculates availability and paired changes. Viral nonsuppression is at least 1,000 copies/mL. Accepted values are 0–5,000 for CD4 count, 0–100 for CD4 percentage, and 0–100 million for viral load. Paired analyses require conception on or after entry and measurements in both phases.

`R/pathway_contacts.R` collapses visits to distinct patient/date contacts. TB encounter flags, chest-radiography entries, or bacteriological-test entries exclude a date from routine contacts. Primary counts exclude the TB event day. Sensitivities vary encounter-code completeness, all-contact inclusion, the 30-day pre-event exclusion, and event-day inclusion.

Uncertainty uses 5,000 woman-clustered bootstrap samples with seed 20260829. Sampled women retain their eligible episodes. `R/manuscript_rng.R` and `src/manuscript_rng.c` reproduce the original PCG64 stream and analysis-specific draw ordering. Generator constants are algorithmic initialization values; no stored participant order or study bootstrap draws are embedded.

`R/pathway_mediation.R` models log(1 + contact dates per 100 observed days) using person-time-weighted linear regression. Poisson outcome models include exposure–mediator interactions. Standardized scenarios use observed person-days and truncate negative transformed-mediator predictions at zero. Direct, indirect, and decomposed total effects use joint woman-clustered covariance and multivariate delta-method intervals. A separately fitted total effect is retained as a separate estimate.

## Tables and figures

`R/table1.R` derives characteristics separately for the cohort and SCCS populations using each design's contributing pregnancies and retained time. Table headings calculate their population sizes from the data.

`R/manuscript_tables.R` creates the manuscript layouts from regenerated aggregate estimates. It preserves the original Wilson–Hilferty approximation for Poisson incidence intervals and applies decimal half-up rounding at presentation. `R/manuscript_documents.R` writes editable Word tables and HTML from the same objects. `R/manuscript_figures.R` draws the illustrative design timeline, pooled-estimate forest plot, participant flow, and postpartum spline.

All observed counts, estimates, intervals, and source-file hashes are computed locally. None are distributed as expected-answer fixtures. Input and intermediate hashes support local provenance and are excluded from code sharing.
