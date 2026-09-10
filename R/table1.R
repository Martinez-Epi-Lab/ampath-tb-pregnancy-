# Baseline characteristics and design-specific pregnancy contributions.
table1_population <- function(cohort, episodes, baseline, retained_rows, ids = cohort$patient_id, design = 'cohort') {
  x<-data.table::copy(data.table::as.data.table(cohort)[patient_id %in% ids])
  b<-data.table::as.data.table(baseline);idx<-match(x$patient_id,b$patient_id)
  x[,`:=`(age_entry=as.numeric(entry-dob)/365.25,art_days_at_entry=as.integer(entry-arvstartTHERAP),
    who=b$who_value[idx],cd4=b$cd4_value[idx],bmi=b$bmi_value[idx],marital=b$marital_value[idx],tpt_before_entry=b$tpt_before_entry[idx])]
  x[,age_category:=c('15–19','20–29','30–39','40–49')[findInterval(age_entry,c(20,30,40))+1L]]
  x[,entry_period:=c('2010–2014','2015–2019','2020–2025')[findInterval(as.integer(format(entry,'%Y')),c(2015,2020))+1L]]
  x[,art_duration_at_entry:=c('<6 months','6–11 months','1–1.99 years','2–4.99 years','≥5 years')[findInterval(art_days_at_entry,c(180,365,730,1825))+1L]]
  x[,`:=`(who_category=data.table::fcase(is.na(who),'Missing',who>=3,'WHO stage 3–4',default='WHO stage 1–2'),
    cd4_category=data.table::fcase(is.na(cd4),'Missing',cd4<200,'<200 cells/µL',cd4<350,'200–349 cells/µL',cd4<500,'350–499 cells/µL',default='≥500 cells/µL'),
    bmi_category=data.table::fcase(is.na(bmi),'Missing',bmi<18.5,'<18.5 kg/m²',bmi<25,'18.5–24.9 kg/m²',default='≥25 kg/m²'),
    partner_category=data.table::fcase(is.na(marital),'Missing',marital==1,'Married/partnered',default='Not married/partnered'))]
  p<-data.table::copy(data.table::as.data.table(episodes)[patient_id %in% x$patient_id & documented_pregnancy_order<=2L])
  ix<-match(p$patient_id,x$patient_id);observation_end<-if(design=='cohort')x$followup_end[ix]else x$administrative_exit[ix]
  p<-p[analysis_start<=observation_end & analysis_end>=x$entry[ix]]
  counts<-p[,.N,by=patient_id];x[,pregnancy_count:=counts$N[match(patient_id,counts$patient_id)]];x[is.na(pregnancy_count),pregnancy_count:=0L]
  x[,`:=`(pregnancy_count_category=as.character(pregnancy_count),multiple_pregnancies=as.integer(pregnancy_count==2L))]
  time<-retained_rows[patient_id %in% ids,.(retained_days=sum(days)),by=patient_id]
  x[,retained_analysis_years:=time$retained_days[match(patient_id,time$patient_id)]/365.25]
  stopifnot(all(x$pregnancy_count<=2L),!anyNA(x$retained_analysis_years))
  x
}

table1_specifications <- function() {
  rows<-list();demographic<-'Demographic and study-entry characteristics';clinical<-'HIV and clinical characteristics at study entry';pregnancy<-'Pregnancy and analysis-specific characteristics'
  add<-function(section,label,column,levels=NULL,digits=1L,binary=FALSE){
    if(is.null(levels))levels<-if(binary)'Yes' else ''
    for(i in seq_along(levels))rows[[length(rows)+1L]]<<-data.table::data.table(section=section,characteristic=if(i==1L)label else '',level=levels[i],column=column,type=if(binary)'binary'else if(levels[i]=='')'continuous'else'category',digits=digits)
  }
  add(demographic,'Age at study entry, years','age_entry');add(demographic,'Age at study entry','age_category',c('15–19','20–29','30–39','40–49'))
  add(demographic,'Calendar period at study entry','entry_period',c('2010–2014','2015–2019','2020–2025'))
  add(demographic,'Time from ART initiation to study entry, days','art_days_at_entry',digits=0L)
  add(demographic,'ART duration at study entry','art_duration_at_entry',c('<6 months','6–11 months','1–1.99 years','2–4.99 years','≥5 years'))
  add(clinical,'WHO stage','who_category',c('WHO stage 1–2','WHO stage 3–4','Missing'))
  add(clinical,'CD4 count among those measured, cells/µL','cd4',digits=0L)
  add(clinical,'CD4 count','cd4_category',c('≥500 cells/µL','350–499 cells/µL','200–349 cells/µL','<200 cells/µL','Missing'))
  add(clinical,'BMI among those measured, kg/m²','bmi',digits=1L)
  add(clinical,'BMI','bmi_category',c('<18.5 kg/m²','18.5–24.9 kg/m²','≥25 kg/m²','Missing'))
  add(clinical,'Partner status','partner_category',c('Married/partnered','Not married/partnered','Missing'))
  add(clinical,'Recorded TPT start date on/before study entry','tpt_before_entry',binary=TRUE)
  add(pregnancy,'Contributing pregnancy episodes among the first two documented','pregnancy_count_category',c('0','1','2'))
  add(pregnancy,'Both first and second documented pregnancies contributed time','multiple_pregnancies',binary=TRUE)
  add(pregnancy,'Retained analysis time after excluding later-pregnancy intervals, years','retained_analysis_years')
  add(pregnancy,'Incident TB event included in the primary analyses','primary_analysis_tb',binary=TRUE)
  out<-data.table::rbindlist(rows);out[,row_number:=seq_len(.N)];out
}

table1_statistic <- function(values,specification) {
  if(specification$type=='continuous'){
    present<-values[!is.na(values)];q<-stats::quantile(present,c(.5,.25,.75),type=7,names=FALSE)
    format<-paste0('%.',specification$digits,'f [%.',specification$digits,'f, %.',specification$digits,'f]')
    return(list(statistic_type='Median [Q1, Q3]',display=sprintf(format,q[1],q[2],q[3]),n_nonmissing=length(present),median=q[1],q1=q[2],q3=q[3]))
  }
  count<-sum(values==if(specification$type=='binary')1L else specification$level,na.rm=TRUE);denominator<-length(values);percent<-count/denominator*100
  list(statistic_type='n (%)',display=sprintf('%s (%.1f%%)',format(count,big.mark=',',trim=TRUE,scientific=FALSE),percent),count=count,denominator=denominator,percent=percent)
}

build_table1_statistics <- function(cohort,sccs) {
  specs<-table1_specifications();groups<-list(cohort=cohort,sccs=sccs)
  data.table::rbindlist(lapply(seq_len(nrow(specs)),function(i){s<-specs[i]
    data.table::rbindlist(lapply(names(groups),function(nm)data.table::as.data.table(c(as.list(s[,.(row_number,section,characteristic,level)]),list(population=nm),table1_statistic(groups[[nm]][[s$column]],s)))),fill=TRUE)
  }),fill=TRUE)
}

table1_formatted <- function(statistics) {
  out<-list();add<-function(label,cohort='',sccs='',type='header',statistic_row=NA_integer_)out[[length(out)+1L]]<<-data.table::data.table(characteristic=label,cohort=cohort,primary_sccs=sccs,row_type=type,statistic_row=statistic_row)
  population_n <- function(group) format(max(statistics[population==group,denominator],na.rm=TRUE),big.mark=',',trim=TRUE)
  add('Characteristic',sprintf('Full cohort (N = %s)',population_n('cohort')),sprintf('Primary SCCS (N = %s)',population_n('sccs')),'column_header');add('','','','blank')
  labels<-c('Median age at study entry, years [IQR]'=1L,'Median time from ART initiation to study entry, days [IQR]'=9L,'Median CD4 count among those measured, cells/µL [IQR]'=18L,'Median BMI among those measured, kg/m² [IQR]'=24L,'Median analysis follow-up time, years [IQR]'=37L)
  specs<-table1_specifications();previous_section<-''
  for(i in seq_len(nrow(specs))){s<-specs[i]
    if(s$section!=previous_section){add(s$section);previous_section<-s$section}
    if(s$type=='category'&&nzchar(s$characteristic))add(if(s$column=='pregnancy_count_category')'Contributing pregnancy episodes'else s$characteristic)
    label<-if(s$type=='category')s$level else s$characteristic
    if(i %in% labels)label<-names(labels)[match(i,labels)]
    if(s$column=='pregnancy_count_category')label<-paste(s$level,if(s$level=='1')'pregnancy'else'pregnancies')
    if(s$column=='primary_analysis_tb')label<-'Incident TB event included in primary analyses'
    add(label,statistics[row_number==i&population=='cohort',display],statistics[row_number==i&population=='sccs',display],'statistic',i)
  }
  add('','','','blank');result<-data.table::rbindlist(out);result[,manuscript_row:=seq_len(.N)-1L];result
}
