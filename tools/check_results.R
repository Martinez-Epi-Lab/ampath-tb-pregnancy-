#!/usr/bin/env Rscript
# Check consistency of locally generated results without stored study values.
args <- commandArgs(FALSE)
self <- gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE)
project <- dirname(dirname(normalizePath(self,mustWork=TRUE)))
for(pkg in c('data.table','jsonlite'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing R package: ',pkg)
report <- file.path(project,'outputs/validation');dir.create(report,recursive=TRUE,showWarnings=FALSE)
checks <- list()
record <- function(name,passed,detail='')checks[[length(checks)+1L]] <<- data.table::data.table(check=name,passed=isTRUE(passed),detail=detail)
figure_names <- c('Figure1_study_design', 'Figure2_primary_IRRs', 'FigureS1_participant_flow', 'FigureS2_postpartum_spline')
for (name in figure_names) for (extension in c('pdf', 'svg', 'png')) {
  path <- file.path(project, 'outputs', 'stage4_figures', paste0(name, '.', extension))
  record(paste0('figure:', name, '.', extension), file.exists(path) && file.info(path)$size > 1000, 'Required nonempty figure artifact')
}
# The generated OOXML has simple, unnested tables and one paragraph per cell.
xml_matches <- function(pattern, text) regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1]]
decode_xml <- function(text) {
  for (entity in c('&lt;', '&gt;', '&quot;', '&apos;', '&amp;')) {
    replacement <- c('&lt;' = '<', '&gt;' = '>', '&quot;' = '"', '&apos;' = "'", '&amp;' = '&')[[entity]]
    text <- gsub(entity, replacement, text, fixed = TRUE)
  }
  text
}
word_cell_text <- function(cell) {
  tokens <- xml_matches('<w:t(?:\\s[^>]*)?>.*?</w:t>|<w:br\\s*/>', cell)
  tokens <- vapply(tokens, function(token) {
    if (grepl('^<w:br', token)) return('\n')
    decode_xml(sub('</w:t>$', '', sub('^<w:t[^>]*>', '', token)))
  }, character(1))
  paste0(tokens, collapse = '')
}
register <- data.table::fread(file.path(project, 'outputs', 'stage4_tables', 'table_register.csv'))
for (name in c('Main_tables', 'Supplementary_tables')) {
  target <- file.path(project, 'outputs', 'stage4_documents', paste0(name, '.docx'))
  expected_count <- if (name == 'Main_tables') 3L else 13L
  tryCatch({
    connection <- unz(target, 'word/document.xml', open = 'r')
    xml <- tryCatch(paste(readLines(connection, warn = FALSE), collapse = ''), finally = close(connection))
    positions <- gregexpr('<w:tbl(?=[ >])', xml, perl = TRUE)[[1]]
    count <- if (positions[1] < 0) 0L else length(positions)
    record(paste0('docx:', name), count == expected_count, sprintf('%d tables; expected %d', count, expected_count))
    tables <- xml_matches('<w:tbl(?:\\s[^>]*)?>.*?</w:tbl>', xml)
    ids <- register$table_id[if (name == 'Main_tables') !grepl('^Table_S', register$table_id) else grepl('^Table_S', register$table_id)]
    if (length(tables) != length(ids)) stop('Word table count does not match table register')
    for (j in seq_along(ids)) {
      expected <- data.table::fread(file.path(project, 'outputs', 'stage4_tables', paste0(ids[j], '.csv')), colClasses = 'character')
      expected_cells <- as.vector(t(rbind(names(expected), as.matrix(expected))))
      cells <- xml_matches('<w:tc(?:\\s[^>]*)?>.*?</w:tc>', tables[j])
      actual_cells <- vapply(cells, word_cell_text, character(1), USE.NAMES = FALSE)
      match <- identical(actual_cells, expected_cells)
      differing <- if (length(actual_cells) == length(expected_cells)) sum(actual_cells != expected_cells, na.rm = TRUE) else NA_integer_
      record(paste0('docx_cells:', ids[j]), match,
        sprintf('%d header/body cells; %s differing cells', length(actual_cells), as.character(differing)))
    }
  }, error = function(e) record(paste0('docx:', name), FALSE, conditionMessage(e)))
}
expected_fits <- c(piece6 = 9L, stage2_secondary = 35L, stage2_spline = 2L, stage3_mediation = 18L)
for (stage in names(expected_fits)) {
  tryCatch({
    diagnostics <- data.table::fread(file.path(project, 'outputs', stage, 'model_diagnostics.csv'))
    parameters <- if ('parameter_count' %in% names(diagnostics)) diagnostics$parameter_count else diagnostics$parameters
    valid <- nrow(diagnostics) == expected_fits[[stage]] && all(!is.na(diagnostics$converged) & diagnostics$converged) &&
      all(!is.na(diagnostics$information_rank) & !is.na(parameters) & diagnostics$information_rank == parameters) &&
      all(is.finite(diagnostics$maximum_absolute_score))
    record(paste0('diagnostics:', stage), valid, sprintf('%d fits; convergence, full information rank and finite scores checked', nrow(diagnostics)))
  }, error = function(e) record(paste0('diagnostics:', stage), FALSE, conditionMessage(e)))
}

flow <- data.table::fread(file.path(project,'outputs/stage4_table1/flow_counts.csv'))
get_count <- function(key) { x<-flow$count[flow$key==key];stopifnot(length(x)==1L);x }
core <- data.table::fread(file.path(project,'outputs/piece6/exposure_estimates.csv'))
for(design in c('cohort','sccs')) {
  index <- which(core$model_id==paste0(design,'_primary'))
  selected <- core[index,]
  record(paste0('event_conservation:',design),sum(selected$events)==get_count('analyzed_events'))
}
record('eligibility_flow',get_count('pregnancy_eligible')-get_count('TB_prior_treatment_excluded')-get_count('TB_broad_evidence_excluded')==get_count('final_cohort'))
record('event_flow',get_count('incident_events')-get_count('excluded_later_pregnancy_events')==get_count('analyzed_events'))
hashes <- data.table::fread(file.path(project,'outputs/stage4_publication/source_checksums.csv'))
record('inputs_unchanged',nrow(hashes)==4L&&all(hashes$unchanged))

result <- data.table::rbindlist(checks)
data.table::fwrite(result,file.path(report,'checks.csv'))
jsonlite::write_json(list(passed=all(result$passed),checks=nrow(result),failed=sum(!result$passed)),file.path(report,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
if(!all(result$passed))stop('Result consistency checks failed; see outputs/validation/checks.csv.')
message('Result consistency checks passed: ',nrow(result),'.')
