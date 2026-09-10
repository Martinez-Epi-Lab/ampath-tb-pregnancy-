#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
self <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE)
project <- dirname(dirname(normalizePath(self,mustWork=TRUE)))
source(file.path(project,'R/share_checks.R'))
check_share(project)
