pathway_phase_order <- function() c('preconception_days_180_to_1','pregnancy','postpartum_days_1_180','postpartum_days_1_30','postpartum_days_31_90','postpartum_days_91_180','postpartum_days_181_365')

build_pathway_episodes <- function(cohort, episodes, patient_order = NULL) {
  c <- data.table::as.data.table(cohort);p <- data.table::copy(data.table::as.data.table(episodes))
  data.table::setorderv(p,c('patient_id','documented_pregnancy_order'))
  p[,next_start:=data.table::shift(start,type='lead')-1L,by=patient_id]
  p<-p[documented_pregnancy_order<=2L]
  idx<-match(p$patient_id,c$patient_id)
  p[,`:=`(entry=c$entry[idx],followup_end=c$followup_end[idx],event_date=c$event_date[idx])]
  p<-p[pmax(analysis_start,entry)<=pmin(analysis_end,followup_end)]
  if(!is.null(patient_order)){p[,patient_rank:=match(patient_id,patient_order)];stopifnot(!anyNA(p$patient_rank));data.table::setorderv(p,c('patient_rank','documented_pregnancy_order'));p[,patient_rank:=NULL]}
  p[,`:=`(episode_number=seq_len(.N),pregnancy_index=documented_pregnancy_order,conceived_after_entry=start>=entry)]
  result<-data.table::rbindlist(lapply(seq_len(nrow(p)),function(i){e<-p[i]
    starts<-c(e$start-180,e$analysis_start,e$end+1,e$end+1,e$end+31,e$end+91,e$end+181)
    ends<-c(e$start-1,e$analysis_end,e$end+180,e$end+30,e$end+90,e$end+180,e$end+365)
    starts<-pmax(starts,e$entry);ends<-pmin(ends,e$followup_end)
    if(!is.na(e$next_start))ends[3:7]<-pmin(ends[3:7],e$next_start)
    out<-data.table::data.table(patient_id=e$patient_id,episode_number=e$episode_number,pregnancy_index=e$pregnancy_index,
      conceived_after_entry=e$conceived_after_entry,pregnancy_start=e$start,pregnancy_end=e$end,
      phase=pathway_phase_order(),phase_start=starts,phase_end=ends)
    out<-out[phase_start<=phase_end];out[,observed_days:=as.integer(phase_end-phase_start)+1L];out
  }))
  result
}

biomarker_histories <- function(visits, ids) {
  v<-data.table::copy(data.table::as.data.table(visits)[patient_id %in% ids])
  valid<-function(x,lo,hi)ifelse(is.finite(x)&x>=lo&x<=hi,x,NA_real_)
  specs<-list(cd4=list(valid(v$cd4,0,5000),v$cd4date),cd4_percentage=list(valid(v$cd4per,0,100),v$cd4date),viral_load=list(valid(v$viralload,0,1e8),v$apptdate),who_stage=list(ifelse(v$WHOstage %in% 1:4,v$WHOstage,NA_real_),v$apptdate))
  bmi<-valid(v$bmi,10,60);bmi[is.na(bmi)]<-valid(v$Weight/(v$Height/100)^2,10,60)[is.na(bmi)]
  specs$bmi<-list(bmi,v$apptdate)
  lapply(specs,function(s){z<-data.table::data.table(patient_id=v$patient_id,date=s[[2]],value=s[[1]],source_row=v$source_row)
    z<-z[!is.na(date)&!is.na(value)];data.table::setorderv(z,c('patient_id','date','source_row'));z[!duplicated(z,by=c('patient_id','date'),fromLast=TRUE)]})
}

biomarker_phase_values <- function(phases,histories) {
  out<-data.table::copy(phases);out[,phase_row:=seq_len(.N)]
  for(nm in names(histories)){
    h<-histories[[nm]]
    picked<-h[out,on=.(patient_id,date>=phase_start,date<=phase_end),nomatch=0,allow.cartesian=TRUE,
      .(phase_row=i.phase_row,measurement_date=x.date,value,source_row)]
    data.table::setorderv(picked,c('phase_row','measurement_date','source_row'))
    summary<-picked[,.(value=value[.N],measurement_date=measurement_date[.N],measurement_dates=.N),by=phase_row]
    ix<-match(out$phase_row,summary$phase_row)
    data.table::set(out,j=nm,value=summary$value[ix]);data.table::set(out,j=paste0(nm,'_date'),value=summary$measurement_date[ix]);data.table::set(out,j=paste0(nm,'_measurement_dates'),value=ifelse(is.na(ix),0L,summary$measurement_dates[ix]))
  }
  out[,phase_row:=NULL];out
}

biomarker_coverage <- function(values) {
  data.table::rbindlist(lapply(pathway_phase_order(),function(ph){x<-values[phase==ph]
    data.table::rbindlist(lapply(c('cd4','cd4_percentage','viral_load','who_stage','bmi'),function(nm){measured<-!is.na(x[[nm]])
      data.table::data.table(phase=ph,candidate=nm,eligible_episodes=nrow(x),eligible_women=data.table::uniqueN(x$patient_id),episodes_measured=sum(measured),coverage_percent=sum(measured)/nrow(x)*100,median_measurement_dates_among_measured=if(any(measured))stats::median(x[[paste0(nm,'_measurement_dates')]][measured]) else NA_real_)}))}))
}

biomarker_pair_frames <- function(values) {
  comparisons<-list(preconception_to_pregnancy=c('preconception_days_180_to_1','pregnancy'),pregnancy_to_postpartum=c('pregnancy','postpartum_days_1_180'),preconception_to_postpartum=c('preconception_days_180_to_1','postpartum_days_1_180'))
  out<-list()
  for(comp in names(comparisons))for(nm in c('cd4','cd4_percentage','viral_load')){
    phases<-comparisons[[comp]];a<-values[conceived_after_entry==TRUE & phase==phases[1],.(patient_id,pregnancy_index,episode_number,before=get(nm))];b<-values[conceived_after_entry==TRUE & phase==phases[2],.(patient_id,pregnancy_index,after=get(nm))]
    z<-merge(a,b,by=c('patient_id','pregnancy_index'));z<-z[!is.na(before)&!is.na(after)];data.table::setorder(z,episode_number)
    z[,`:=`(comparison=comp,candidate=nm,change=after-before)]
    if(nm=='viral_load')z[,`:=`(before_binary=as.integer(before>=1000),after_binary=as.integer(after>=1000),change=100*(as.integer(after>=1000)-as.integer(before>=1000)))]
    out[[paste(comp,nm,sep=':')]]<-z
  }
  out
}

biomarker_bootstrap <- function(pairs,rng,replicates=5000L) {
  codes<-match(pairs$patient_id,unique(pairs$patient_id));n<-max(codes);values<-pairs$change
  ord<-order(values,seq_along(values));sv<-values[ord];sc<-codes[ord]
  estimates<-numeric(replicates)
  for(first in seq.int(1L,replicates,by=100L)){
    number<-min(100L,replicates-first+1L);counts<-manuscript_rng_counts(rng,n,number,'integers')
    for(j in seq_len(number)){
      if(pairs$candidate[1]=='viral_load'){w<-counts[codes,j];estimates[first+j-1L]<-sum(w*values)/sum(w)}
      else {w<-counts[sc,j];cum<-cumsum(w);total<-sum(w);left<-(total-1L)%/%2L+1L;right<-total%/%2L+1L;estimates[first+j-1L]<-(sv[which(cum>=left)[1]]+sv[which(cum>=right)[1]])/2}
    }
  }
  list(ci=stats::quantile(estimates,c(.025,.975),type=7,names=FALSE),replicates=estimates)
}

summarize_biomarker_pairs <- function(frames,rng,replicates=5000L) {
  draws<-list();summaries<-list()
  for(nm in names(frames)){
    p<-frames[[nm]];if(!nrow(p))next
    boot<-biomarker_bootstrap(p,rng,replicates);draws[[nm]]<-boot$replicates
    r<-list(comparison=p$comparison[1],candidate=p$candidate[1],paired_episodes=nrow(p),paired_women=data.table::uniqueN(p$patient_id),change_type=if(p$candidate[1]=='viral_load')'paired prevalence difference, percentage points' else 'median paired within-episode change',change_estimate=if(p$candidate[1]=='viral_load')mean(p$change) else stats::median(p$change),change_ci_lower=boot$ci[1],change_ci_upper=boot$ci[2])
    if(p$candidate[1]=='viral_load')r<-c(r,list(earlier_n=sum(p$before_binary),earlier_percent=mean(p$before_binary)*100,later_n=sum(p$after_binary),later_percent=mean(p$after_binary)*100,changed_to_nonsuppressed_n=sum(p$before_binary==0&p$after_binary==1),changed_to_suppressed_n=sum(p$before_binary==1&p$after_binary==0)))
    else {q<-function(x)stats::quantile(x,c(.5,.25,.75),type=7,names=FALSE);r<-c(r,as.list(stats::setNames(c(q(p$before),q(p$after),q(p$change)[2:3]),c('earlier_median','earlier_q1','earlier_q3','later_median','later_q1','later_q3','change_q1','change_q3'))))}
    summaries[[nm]]<-data.table::as.data.table(r)
  }
  list(summary=data.table::rbindlist(summaries,fill=TRUE),bootstrap_replicates=draws)
}
