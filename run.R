#!/usr/bin/env Rscript
# Run the manuscript analysis from source extracts to tables and figures.
args <- commandArgs(FALSE)
self <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE)
project <- dirname(normalizePath(self,mustWork=TRUE))
source(file.path(project,'R/pipeline.R'))
steps <- pipeline_steps()
options <- commandArgs(TRUE)
if(identical(options,'--list')) { print(steps[,c('number','name','description')],row.names=FALSE); quit(status=0) }
if(any(!grepl('^--(from|through)=[0-9]+$',options))) stop('Usage: Rscript --vanilla run.R [--from=1] [--through=18] or --list')
get_number <- function(key,default) {
  value <- options[startsWith(options,paste0('--',key,'='))]
  if(length(value)>1L)stop('Repeated option: ',key)
  if(length(value))as.integer(sub('.*=','',value))else default
}
first <- get_number('from',1L); last <- get_number('through',nrow(steps))
if(anyNA(c(first,last))||first<1L||last>nrow(steps)||first>last)stop('Invalid module range.')
for(pkg in c('haven','data.table','jsonlite','digest','survival'))
  if(!requireNamespace(pkg,quietly=TRUE))stop('Missing R package: ',pkg)
if(!nzchar(Sys.which('zip')))stop('The zip command is required.')
if(!any(nzchar(Sys.which(c('clang','cc')))))stop('A C compiler is required.')
config <- file.path(project,'config/paths.local.R')
if(!file.exists(config))stop('Copy config/paths.example.R to config/paths.local.R and set all four input paths.')
cfg <- new.env(parent=baseenv());sys.source(config,envir=cfg)
paths <- c(unlist(cfg$ampath_paths,use.names=FALSE),cfg$ampath_outcome_path)
if(length(paths)!=4L||anyNA(paths)||any(!file.exists(paths)))stop('All four configured extracts must exist.')
for(path in normalizePath(paths,mustWork=TRUE))
  if(startsWith(path,paste0(project,'/')))stop('Keep the source extracts outside the code project.')
run_script <- function(path) {
  status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(file.path(project,path))))
  if(status!=0L)stop('Analysis stopped at ',path,'.')
}
for(i in seq.int(first,last)) {
  message(sprintf('[%02d/%02d] %s',i,nrow(steps),steps$description[i]))
  run_script(steps$script[i])
}
if(last==nrow(steps))run_script('tools/check_results.R')
message('Requested analysis modules completed. Results are local under outputs/.')
