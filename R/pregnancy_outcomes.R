# Outcome documentation is linked to original pregnancy ends before adjudication.
outcome_nearest_code <- function(x, codes) {
  z <- x[pregnancyoutcome %in% codes]
  if (!nrow(z)) return(NA_integer_)
  winner <- unique(z[abs(days_from_pregnancy_end) == min(abs(days_from_pregnancy_end)), pregnancyoutcome])
  if (length(winner) == 1L) winner else NA_integer_
}

adjudicate_pregnancy_outcome <- function(records, max_after = 90L) {
  x <- records[days_from_pregnancy_end >= -30 & days_from_pregnancy_end <= max_after]
  prefix <- if (max_after == 90L) "primary90" else "extended180"
  result <- function(count, codes, vital, gestation, neonatal, summary, adverse, conflict, confidence)
    stats::setNames(list(count, length(codes), paste(sort(codes), collapse = ","), vital,
      gestation, neonatal, summary, adverse, conflict, confidence), paste0(prefix, "_",
      c("record_count", "distinct_code_count", "codes", "vital_outcome", "gestational_category",
        "neonatal_outcome", "summary_category", "any_adverse_recorded", "conflict_status", "adjudication_confidence")))
  if (!nrow(x)) return(result(0L, integer(), "Unknown", "Unknown", "Unknown", "Unknown outcome", NA,
    "No linked outcome record", "Unknown"))
  codes <- unique(x$pregnancyoutcome)
  direct <- intersect(codes, c(1L,4L,5L)); implied <- intersect(codes,8L); loss <- intersect(codes,c(2L,3L))
  notes <- character(); resolved <- unresolved <- FALSE
  note <- function(text, unres = FALSE) { notes <<- c(notes,text); if (unres) unresolved <<- TRUE else resolved <<- TRUE }
  if (length(loss) && length(direct)) {
    winner <- outcome_nearest_code(x,c(loss,direct))
    if (is.na(winner)) { vital <- "Unresolved live-birth/pregnancy-loss conflict"; note("live_vs_loss_unresolved_exact_timing_tie",TRUE) }
    else { vital <- if (winner==2L) "Miscarriage" else if (winner==3L) "Stillbirth" else "Live birth"; note("live_vs_loss_resolved_by_nearest_explicit_code") }
  } else if (length(loss)) {
    if (all(c(2L,3L) %in% loss)) {
      winner <- outcome_nearest_code(x,c(2L,3L))
      if (is.na(winner)) { vital <- "Pregnancy loss; type unresolved"; note("miscarriage_vs_stillbirth_unresolved_exact_timing_tie",TRUE) }
      else { vital <- if (winner==2L) "Miscarriage" else "Stillbirth"; note("miscarriage_vs_stillbirth_resolved_by_nearest_code") }
    } else {
      vital <- if (2L %in% loss) "Miscarriage" else "Stillbirth"
      if (length(implied)) note("loss_vs_term_infant_resolved_by_clinical_specificity")
    }
  } else if (length(direct) || length(implied)) vital <- "Live birth"
  else vital <- if (7L %in% codes) "Delivery; vital status not recorded" else "Unknown"
  if (vital=="Miscarriage") gestation <- "Not classified after miscarriage"
  else if (all(c(7L,8L) %in% codes)) {
    winner <- outcome_nearest_code(x,c(7L,8L))
    if (is.na(winner)) { gestation <- "Unresolved preterm/term conflict";note("preterm_vs_term_unresolved_exact_timing_tie",TRUE) }
    else {gestation <- if(winner==7L) "Preterm" else "Term";note("preterm_vs_term_resolved_by_nearest_code")}
  } else gestation <- if(7L %in% codes) "Preterm" else if(8L %in% codes) "Term" else "Unknown"
  if (all(c(4L,5L) %in% codes)) {
    winner <- outcome_nearest_code(x,c(4L,5L))
    if(is.na(winner)) {neonatal<-"Neonatal death; timing category unresolved";note("early_vs_late_neonatal_death_unresolved_exact_timing_tie",TRUE)}
    else {neonatal<-if(winner==4L) "Neonatal death within 7 days" else "Neonatal death at 7-29 days";note("early_vs_late_neonatal_death_resolved_by_nearest_code")}
  } else neonatal <- if(4L %in% codes) "Neonatal death within 7 days" else if(5L %in% codes) "Neonatal death at 7-29 days" else if(vital=="Live birth") "No neonatal death recorded" else "Unknown"
  if(vital %in% c("Miscarriage","Stillbirth") || startsWith(vital,"Unresolved") || vital=="Pregnancy loss; type unresolved") summary<-vital
  else if(vital=="Live birth") summary<-if(startsWith(neonatal,"Neonatal death")) "Live birth with neonatal death" else if(gestation=="Preterm") "Live birth, preterm" else if(gestation=="Term") "Live birth, term" else if(startsWith(gestation,"Unresolved")) "Live birth; preterm/term status unresolved" else "Live birth; gestational category not recorded"
  else summary<-if(gestation=="Preterm") "Preterm delivery; vital status not recorded" else "Outcome information insufficient"
  adverse<-if(vital %in% c("Miscarriage","Stillbirth") || startsWith(neonatal,"Neonatal death") || gestation=="Preterm") TRUE else if(unresolved) NA else FALSE
  tier<-min(ifelse(x$days_from_pregnancy_end<=30,1L,ifelse(x$days_from_pregnancy_end<=90,2L,3L)))
  confidence<-if(unresolved) "Unresolved" else if(resolved || tier==3L) "Low" else if(tier==2L) "Moderate" else "High"
  result(nrow(x),codes,vital,gestation,neonatal,summary,adverse,if(length(notes)) paste(notes,collapse="; ") else "No contradictory codes",confidence)
}

link_pregnancy_outcome_records <- function(records, episodes) {
  e <- data.table::copy(data.table::as.data.table(episodes)[documented_pregnancy_order<=2,
    .(patient_id, documented_pregnancy_order, original_end)])
  stopifnot(!anyDuplicated(e[,.(patient_id,documented_pregnancy_order)]))
  candidates <- merge(records,e,by="patient_id",allow.cartesian=TRUE)
  candidates[,days_from_pregnancy_end:=as.integer(outcome_record_date-original_end)]
  candidates <- candidates[!is.na(days_from_pregnancy_end) & days_from_pregnancy_end>=-30 & days_from_pregnancy_end<=180]
  candidates[,`:=`(absolute_days_from_end=abs(days_from_pregnancy_end),before_tie_break=as.integer(days_from_pregnancy_end<0))]
  data.table::setorderv(candidates,c("outcome_record_id","absolute_days_from_end","before_tie_break","documented_pregnancy_order"))
  assignments <- candidates[!duplicated(outcome_record_id)]
  groups <- split(assignments,paste(assignments$patient_id,assignments$documented_pregnancy_order,sep="|"))
  empty <- assignments[0]
  adjudicated <- data.table::rbindlist(lapply(seq_len(nrow(e)),function(i) {
    x<-groups[[paste(e$patient_id[i],e$documented_pregnancy_order[i],sep="|")]];if(is.null(x))x<-empty
    data.table::as.data.table(c(as.list(e[i]),adjudicate_pregnancy_outcome(x,90L),adjudicate_pregnancy_outcome(x,180L)))
  }),fill=TRUE)
  adjudicated[,recorded_live_birth:=primary90_vital_outcome=="Live birth"]
  # Repeated codes remain in record counts; collapsed evidence is a separate audit.
  collapsed <- assignments[,.(source_record_count=data.table::uniqueN(outcome_record_id),
    first_record_date=min(outcome_record_date),last_record_date=max(outcome_record_date)),
    by=.(patient_id,documented_pregnancy_order,pregnancyoutcome)]
  list(outcomes=adjudicated,assignments=assignments,collapsed_evidence=collapsed)
}

derive_pregnancy_outcomes <- function(episodes, visits, outcome_path) {
  path <- normalizePath(outcome_path,mustWork=TRUE)
  before <- digest::digest(path,algo="sha256",file=TRUE)
  raw <- data.table::fread(path,select=c("patient_id","apptdate","pregnancyoutcome"),colClasses="character",na.strings=NULL,strip.white=FALSE)
  v <- data.table::copy(data.table::as.data.table(visits)[,.(source_row,patient_id,source_apptdate)])
  data.table::setorder(v,source_row)
  stopifnot(nrow(raw)==nrow(v),identical(v$source_row,seq_len(nrow(v))))
  ids <- canonical_id(raw$patient_id)
  date_tokens <- blank_to_na(raw$apptdate)
  if(!identical(ids,v$patient_id) || !identical(date_tokens,blank_to_na(v$source_apptdate))) stop("Outcome extract patient/date row keys do not exactly match visits.")
  codes <- parse_numeric_field(raw$pregnancyoutcome)
  if(any(!is.na(codes) & (codes!=floor(codes) | !codes %in% 1:8))) stop("Unrecognized pregnancy outcome code.")
  keep <- !is.na(codes) & ids %in% episodes$patient_id
  records <- data.table::data.table(patient_id=ids[keep],outcome_record_id=seq_len(sum(keep)),
    visit_source_row=which(keep),outcome_record_date=parse_ampath_date(raw$apptdate[keep]),pregnancyoutcome=as.integer(codes[keep]))
  result <- link_pregnancy_outcome_records(records,episodes)
  after <- digest::digest(path,algo="sha256",file=TRUE);stopifnot(identical(before,after))
  result$aggregateQC <- list(source_rows=nrow(raw),source_sha256=before,source_unchanged=TRUE,
    patient_key_mismatches=0L,date_key_mismatches=0L,all_outcome_rows=sum(!is.na(codes)),
    all_outcome_women=data.table::uniqueN(ids[!is.na(codes)]),analytic_outcome_rows=nrow(records),
    analytic_outcome_women=data.table::uniqueN(records$patient_id),assigned_records=nrow(result$assignments),
    first_two_episodes=nrow(result$outcomes),primary_recorded_episodes=sum(result$outcomes$primary90_record_count>0),
    recorded_live_birth_episodes=sum(result$outcomes$recorded_live_birth))
  result
}
