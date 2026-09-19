#!/usr/bin/env Rscript
args<-commandArgs(FALSE);this_file<-gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE)
project <- dirname(dirname(normalizePath(this_file,mustWork=TRUE)));suppressPackageStartupMessages(library(data.table))
for(f in c('stage_io.R','stage2_io.R','manuscript_rng.R','pathway_biomarkers.R'))source(file.path(project,'R',f))
ctx<-start_stage(project,'stage3_biomarkers');check_core_versions(project);load_manuscript_rng(project)
inputs<-outputs<-list()
read_input<-function(stage,file,key){x<-stage2_input(ctx,stage,file);inputs[[key]]<<-x$fingerprint;x$data}
c<-read_input('piece4','cohort.rds','cohort');p<-read_input('piece4','episodes.rds','episodes');a<-read_input('piece3','audit.rds','pregnancy_audit');allp<-read_input('piece3','episodes.rds','all_episodes');v<-read_input('piece1','visits.rds','visits')
patient_order<-python_patient_order(a,allp)
message('Assigning biomarker measurements to pregnancy phases.')
phases<-build_pathway_episodes(c,p,patient_order);histories<-biomarker_histories(v,c$patient_id);values<-biomarker_phase_values(phases,histories)
# Manuscript summaries use preconception, pregnancy, and postpartum days 1–180.
values<-values[phase %in% pathway_phase_order()[1:3]]
pairs<-biomarker_pair_frames(values)
outputs$phase_values.rds<-save_stage_rds(values,ctx$work,'phase_values.rds')
outputs$paired_episodes.rds<-save_stage_rds(pairs,ctx$work,'paired_episodes.rds')
coverage<-biomarker_coverage(values)[phase %in% pathway_phase_order()[1:3]&candidate %in% c('cd4','cd4_percentage','viral_load')]
fwrite(coverage,file.path(ctx$report,'biomarker_phase_coverage.csv'))
message('Running 5,000 bootstrap replicates, resampling women, for nine comparisons.')
paired<-summarize_biomarker_pairs(pairs,new_manuscript_rng(),5000L)
outputs$bootstrap_replicates.rds<-save_stage_rds(paired$bootstrap_replicates,ctx$work,'bootstrap_replicates.rds')
fwrite(paired$summary,file.path(ctx$report,'biomarker_paired_changes.csv'))
summary<-list(phase_rows=nrow(values),episodes=uniqueN(values$episode_number),women=uniqueN(values$patient_id),paired_comparisons=nrow(paired$summary),bootstrap_replicates=5000L)
finish_stage(ctx,inputs,outputs,list(seed=20260829L,generator='NumPy-compatible PCG64 with bounded-integer sampling',bootstrap='Resample women, retaining all paired episodes; percentile CI, linear quantiles',ordering='CPython3.12 set(defaultdict) order reconstructed from source',phase_value='Latest valid measurement within phase',TB_censoring='Inclusive first incident TB event',phase_names=pathway_phase_order()[1:3]),summary,c('R/stage2_io.R','R/pathway_biomarkers.R','R/manuscript_rng.R','src/manuscript_rng.c'),list(prior_TB_exclusions_applied=TRUE,TB_events_ascertained=TRUE,first_two_restriction_applied=TRUE))
jsonlite::write_json(summary,file.path(ctx$report,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
message('Biomarker analyses complete.')
