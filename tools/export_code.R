#!/usr/bin/env Rscript
# Export only the reviewed code and documentation; no data or Git history.
args <- commandArgs(FALSE)
self <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE)
project <- dirname(dirname(normalizePath(self,mustWork=TRUE)))
source(file.path(project,'R/share_checks.R'))
files <- check_share(project)
if(!nzchar(Sys.which('zip')))stop('The zip command is required.')
release <- file.path(project,'release');dir.create(release,showWarnings=FALSE)
scratch <- tempfile('code_export_');dir.create(scratch)
on.exit_cleanup <- function()unlink(scratch,recursive=TRUE)
tryCatch({
  folder <- file.path(scratch,'ampath-tb-pregnancy');dir.create(folder)
  for(file in files) {
    target <- file.path(folder,file);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(file.path(project,file),target))stop('Could not copy ',file)
  }
  check_share(folder)
  Sys.setFileTime(file.path(folder,files),as.POSIXct('2000-01-01',tz='UTC'))
  archive <- file.path(release,'ampath-tb-pregnancy-code-only.zip')
  if(file.exists(archive))unlink(archive)
  old <- getwd();setwd(scratch)
  status <- tryCatch(utils::zip(archive,file.path('ampath-tb-pregnancy',sort(files)),flags='-q -X'),finally=setwd(old))
  if(status!=0L)stop('Code archive creation failed.')
  inventory <- utils::unzip(archive,list=TRUE)$Name
  if(!identical(sort(inventory),sort(file.path('ampath-tb-pregnancy',files))))stop('Archive does not match approved files.')
  message('Created release/ampath-tb-pregnancy-code-only.zip with ',length(files),' files.')
},finally=on.exit_cleanup())
