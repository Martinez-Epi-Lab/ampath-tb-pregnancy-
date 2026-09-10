#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
self <- gsub("~+~", " ", sub("^--file=", "", args[grepl("^--file=", args)][1]), fixed = TRUE)
project <- dirname(dirname(normalizePath(self, mustWork = TRUE)))
for (f in c("stage_io.R", "stage2_io.R", "manuscript_documents.R")) source(file.path(project, "R", f))
ctx <- start_stage(project, "stage4_documents")
input <- stage2_input(ctx, "stage4_tables", "tables.rds"); tables <- input$data
main <- tables[!startsWith(names(tables), "Table_S")]; supplement <- tables[startsWith(names(tables), "Table_S")]
write_manuscript_docx(main, file.path(ctx$report, "Main_tables.docx"), file.path(ctx$work, "main_ooxml"))
write_manuscript_docx(supplement, file.path(ctx$report, "Supplementary_tables.docx"), file.path(ctx$work, "supplement_ooxml"))
write_manuscript_html(tables, file.path(ctx$report, "Manuscript_tables.html"))
paths <- file.path(ctx$report, c("Main_tables.docx", "Supplementary_tables.docx", "Manuscript_tables.html"))
hashes <- setNames(lapply(paths, digest::digest, file=TRUE, algo="sha256"), basename(paths))
finish_stage(ctx, list(tables=input$fingerprint), hashes, list(editable_Word_tables=TRUE, source="Regenerated R table objects"),
  list(main_tables=length(main),supplementary_panels=length(supplement)), "R/manuscript_documents.R")
message("Wrote editable main and supplementary Word tables and HTML view.")
