#!/usr/bin/env Rscript
args<-commandArgs(FALSE);this_file<-gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE)
project <- dirname(dirname(normalizePath(this_file,mustWork=TRUE)));suppressPackageStartupMessages(library(data.table))
for(f in c('stage_io.R','manuscript_figures.R'))source(file.path(project,'R',f))
ctx<-start_stage(project,'stage4_figures');inputs<-outputs<-list()
read_aggregate<-function(stage,file,key){path<-file.path(project,'outputs',stage,file);inputs[[key]]<<-list(path=path,sha256=digest::digest(path,file=TRUE,algo='sha256'));fread(path)}
coefficients<-read_aggregate('piece6','model_coefficients.csv','primary_coefficients')
flow<-read_aggregate('stage4_table1','flow_counts.csv','flow_counts')
curves<-read_aggregate('stage2_spline','postpartum_spline_curve.csv','spline_curves');tests<-read_aggregate('stage2_spline','postpartum_spline_tests.csv','spline_tests')
figures<-list(Figure1_study_design=list(width=9,height=6.1,draw=function()figure_study_design()),Figure2_primary_IRRs=list(width=8.2,height=4.8,draw=function()figure_primary_forest(coefficients)),FigureS1_participant_flow=list(width=8.4,height=9.2,draw=function()figure_participant_flow(flow)),FigureS2_postpartum_spline=list(width=10.5,height=5.2,draw=function()figure_postpartum_spline(curves,tests)))
for(name in names(figures)){spec<-figures[[name]];message('Drawing ',name,'.');paths<-save_manuscript_figure(ctx$report,name,spec$width,spec$height,spec$draw);for(path in paths)outputs[[basename(path)]]<-digest::digest(path,file=TRUE,algo='sha256')}
manifest<-finish_stage(ctx,inputs,outputs,list(formats=c('PDF','SVG','PNG'),png_dpi=300L,figure1='Illustrative timeline; no participant records',forest='Canonical age/calendar-adjusted primary estimates',flow='Counts recomputed from core audits',spline='Daily fitted curves and95%intervals from validated stage2 outputs'),list(figures=length(figures),files=length(outputs)),c('R/manuscript_figures.R'),list(prior_TB_exclusions_applied=TRUE,TB_events_ascertained=TRUE,first_two_restriction_applied=TRUE,models_fitted=TRUE))
fwrite(data.table(figure=names(figures),pdf=paste0(names(figures),'.pdf'),svg=paste0(names(figures),'.svg'),png=paste0(names(figures),'.png')),file.path(ctx$report,'figure_inventory.csv'))
message('Four manuscript figures generated.')
