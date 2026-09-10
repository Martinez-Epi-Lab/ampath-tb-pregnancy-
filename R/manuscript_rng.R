# The manuscript used NumPy PCG64, with seed 20260829.
load_manuscript_rng <- function(project) {
  source<-file.path(project,'src','manuscript_rng.c');dir<-file.path(project,'work','native');dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  library<-file.path(dir,paste0('manuscript_rng',.Platform$dynlib.ext))
  hash<-digest::digest(source,algo='sha256',file=TRUE);stamp<-paste0(library,'.sha256')
  if(!file.exists(library)||!file.exists(stamp)||readLines(stamp,warn=FALSE)!=hash){
    cc<-Sys.which('clang');if(!nzchar(cc))cc<-Sys.which('cc');if(!nzchar(cc))stop('A C compiler is required for the manuscript RNG.')
    args<-c('-O2','-fPIC','-shared',if(Sys.info()[['sysname']]=='Darwin')c('-undefined','dynamic_lookup'),paste0('-I',shQuote(R.home('include'))),shQuote(source),'-o',shQuote(library))
    status<-system2(cc,args);if(status!=0)stop('Could not build manuscript RNG.');writeLines(hash,stamp)
  }
  if(!is.loaded('ampath_rng_new'))dyn.load(library)
  invisible(library)
}
new_manuscript_rng <- function(seed=20260829L) .Call('ampath_rng_new',as.integer(seed))
manuscript_rng_counts <- function(rng,n,replicates,method=c('integers','multinomial')) .Call('ampath_rng_counts',rng,as.integer(n),as.integer(replicates),match.arg(method)=='multinomial')
python_patient_order <- function(pregnancy_audit,clean_episodes) {
  p<-data.table::copy(data.table::as.data.table(pregnancy_audit)[pregnancy_status %in% c('retained','collapsed_after_overlap','empty_after_observation_clipping')])
  data.table::setorder(p,pregnancy_source_row)
  ids<-unique(p$patient_id);ids<-ids[ids %in% clean_episodes$patient_id]
  ids[.Call('ampath_python_set_order',as.numeric(ids))]
}
