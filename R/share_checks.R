# Check files listed in share-files.txt before export.
share_files <- function(project) {
  files <- readLines(file.path(project,'share-files.txt'),warn=FALSE)
  files <- files[nzchar(files)]
  if(!length(files)||anyDuplicated(files)||any(grepl('(^/|\\\\|(^|/)\\.\\.(/|$))',files)))stop('Invalid sharing file list.')
  files
}
check_share <- function(project) {
  project <- normalizePath(project,mustWork=TRUE)
  files <- share_files(project)
  failures <- character()
  fail <- function(file,reason)failures <<- c(failures,paste(file,reason,sep=': '))
  approved_extension <- '[.](R|c|md|txt|yml|Rproj)$|^([.]gitignore|[.]gitattributes)$'
  patterns <- c(personal_path='/(?:Users|home)/[^/[:space:]"<>]+|/[Pp]rivate/(?:tmp|var)/',
    cloud_account='Google[D]rive-|Cloud[S]torage/|shortcut-targets-by-[i]d',
    email='[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}',
    private_key='-----BEGIN [A-Z ]*PRIVATE KEY-----',
    credential='(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9_-]{24,}|AKIA[A-Z0-9]{16})')
  for(file in files) {
    path <- file.path(project,file)
    if(!grepl(approved_extension,file)){fail(file,'unapproved file type');next}
    if(!file.exists(path)||dir.exists(path)){fail(file,'missing file');next}
    ancestors <- path
    while(dirname(tail(ancestors,1))!=project)ancestors<-c(ancestors,dirname(tail(ancestors,1)))
    if(any(nzchar(Sys.readlink(ancestors)),na.rm=TRUE)){fail(file,'symbolic link');next}
    if(file.info(path)$size>10*1024^2){fail(file,'unexpectedly large source file');next}
    bytes <- readBin(path,'raw',n=file.info(path)$size)
    if(any(bytes==as.raw(0))){fail(file,'binary content');next}
    text <- rawToChar(bytes)
    if(is.na(iconv(text,from='UTF-8',to='UTF-8'))){fail(file,'invalid UTF-8');next}
    for(name in names(patterns))if(grepl(patterns[[name]],text,perl=TRUE))fail(file,paste('possible',name))
  }
  all_files <- list.files(project,recursive=TRUE,all.files=TRUE,include.dirs=FALSE,no..=TRUE)
  private_local <- grepl('^(work|outputs|release|data|raw|[.]git|[.]runtime)/|^config/paths[.]local[.]R$',all_files)
  unexpected <- setdiff(all_files[!private_local],files)
  for(file in unexpected)fail(file,'not approved for sharing')
  if(nzchar(Sys.which('git'))) {
    top <- suppressWarnings(system2('git',c('-C',shQuote(project),'rev-parse','--show-toplevel'),stdout=TRUE,stderr=FALSE))
    if(length(top)==1L&&identical(normalizePath(top,mustWork=FALSE),project)) {
      tracked <- system2('git',c('-C',shQuote(project),'ls-files'),stdout=TRUE)
      for(file in setdiff(tracked,files))fail(file,'unapproved tracked file')
      for(file in intersect(tracked,files)) {
        staged <- system2('git',c('-C',shQuote(project),'show',shQuote(paste0(':',file))),stdout=TRUE)
        text <- paste(staged,collapse='\n')
        for(name in names(patterns))if(grepl(patterns[[name]],text,perl=TRUE))fail(file,paste('possible staged',name))
      }
    }
  }
  if(length(failures))stop(paste(c('Sharing check failed:',failures),collapse='\n'))
  message('Sharing check passed: ',length(files),' approved text files.')
  invisible(files)
}
