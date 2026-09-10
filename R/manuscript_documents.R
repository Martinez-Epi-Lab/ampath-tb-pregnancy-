# Editable Word tables and an HTML view from the same result objects.
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE)
}
word_paragraph <- function(text, style = "Normal", bold = FALSE, size = 20L, align = "left", keep = FALSE) {
  pieces <- strsplit(as.character(text), "\n", fixed = TRUE)[[1]]
  if (!length(pieces)) pieces <- ""
  content <- paste0('<w:t xml:space="preserve">', xml_escape(pieces), '</w:t>', collapse = '<w:br/>')
  paste0('<w:p><w:pPr><w:pStyle w:val="', style, '"/><w:jc w:val="', align, '"/>',
    if (keep) '<w:keepNext/>' else '', '<w:spacing w:after="35" w:line="240" w:lineRule="auto"/></w:pPr><w:r><w:rPr>',
    if (bold) '<w:b/>' else '', '<w:sz w:val="', size, '"/><w:szCs w:val="', size, '"/></w:rPr>', content, '</w:r></w:p>')
}
word_column_weights <- function(id, n) {
  specs <- list(Table_1 = c(.60,.20,.20), Table_2 = c(.40,.13,.16,.31), Table_3 = c(.34,rep(.165,4)),
    Table_S1 = c(.16,.44,.20,.20), Table_S2 = c(.20,rep(.10,8)), Table_S3 = c(.07,.20,.10,.07,.10,.23,.23),
    Table_S4 = c(.10,.25,.325,.325), Table_S5 = c(.08,.14,.10,.22,.14,.10,.22), Table_S6 = c(.15,.13,.13,.19,.20,.20),
    Table_S7A = c(.19,.17,.14,.15,.15,.20), Table_S7B = c(.18,.13,.17,.17,.15,.20),
    Table_S8 = c(.23,.13,.17,.20,.27), Table_S9 = c(.24,.16,.19,.19,.22),
    Table_S10A = c(.15,.13,rep(.18,4)), Table_S10B = c(.46,.27,.27), Table_S11 = c(.09,.12,.285,.285,.22))
  weights <- specs[[id]]; if (is.null(weights)) weights <- c(2, rep(1,n-1))
  weights / sum(weights)
}
word_table <- function(tab) {
  d <- tab$data; total <- if (tab$landscape) 14112L else 10512L
  widths <- as.integer(total * word_column_weights(tab$id, ncol(d))); widths[length(widths)] <- total - sum(head(widths,-1))
  grid <- paste0('<w:gridCol w:w="', widths, '"/>', collapse = '')
  borders <- '<w:tblBorders><w:top w:val="single" w:sz="10" w:color="000000"/><w:bottom w:val="single" w:sz="10" w:color="000000"/><w:left w:val="nil"/><w:right w:val="nil"/><w:insideH w:val="nil"/><w:insideV w:val="nil"/></w:tblBorders>'
  rows <- c(list(names(d)), lapply(seq_len(nrow(d)), function(i) as.character(d[i,])))
  body <- vapply(seq_along(rows), function(i) {
    header <- i == 1L; group <- (i - 1L) %in% tab$groups
    cells <- vapply(seq_along(widths), function(j) {
      shade <- if (header) 'F1F1F1' else if (group) 'F7F7F7' else 'FFFFFF'
      bottom <- if (header) '<w:tcBorders><w:bottom w:val="single" w:sz="6" w:color="000000"/></w:tcBorders>' else ''
      left <- j == 1L || grepl('^(Adjustment|Comparison|Measure|Exposure|Period|Contact definition)', names(d)[j])
      paste0('<w:tc><w:tcPr><w:tcW w:w="', widths[j], '" w:type="dxa"/><w:shd w:fill="', shade, '"/><w:vAlign w:val="center"/>', bottom,
        '</w:tcPr>', word_paragraph(rows[[i]][j], bold = header || group, size = if (tab$id == 'Table_1') 19L else 20L,
          align = if (left) 'left' else 'center'), '</w:tc>')
    }, character(1))
    paste0('<w:tr><w:trPr><w:cantSplit/>', if (header) '<w:tblHeader/>' else '', '</w:trPr>', paste0(cells, collapse=''), '</w:tr>')
  }, character(1))
  paste0('<w:tbl><w:tblPr><w:tblW w:w="',total,'" w:type="dxa"/><w:tblLayout w:type="fixed"/>', borders,
    '<w:tblCellMar><w:top w:w="48" w:type="dxa"/><w:bottom w:w="48" w:type="dxa"/><w:left w:w="80" w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid>', grid, '</w:tblGrid>', paste0(body,collapse=''), '</w:tbl>')
}
word_section <- function(landscape, paragraph = FALSE) {
  w <- if (landscape) 15840 else 12240; h <- if (landscape) 12240 else 15840
  x <- paste0('<w:sectPr><w:footerReference w:type="default" r:id="rIdFooter"/><w:type w:val="nextPage"/><w:pgSz w:w="',w,'" w:h="',h,'"',
    if (landscape) ' w:orient="landscape"' else '', '/><w:pgMar w:top="864" w:right="864" w:bottom="864" w:left="864" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>')
  if (paragraph) paste0('<w:p><w:pPr>',x,'</w:pPr></w:p>') else x
}
write_manuscript_docx <- function(tables, output, scratch) {
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  for (dir in c('_rels','word/_rels','docProps')) dir.create(file.path(scratch,dir),recursive=TRUE,showWarnings=FALSE)
  write_xml <- function(path,text) writeLines(paste0('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',text),file.path(scratch,path),useBytes=TRUE)
  ids <- names(tables); body <- character()
  for (i in seq_along(ids)) {
    tab <- tables[[i]]
    body <- c(body, word_paragraph(tab$title, 'Title', TRUE, 28L, keep=TRUE), word_table(tab), word_paragraph(tab$notes, 'Note', FALSE, 18L))
    same_next <- i < length(ids) && sub('[AB]$','',ids[i]) == sub('[AB]$','',ids[i+1L])
    if (i < length(ids) && !same_next) body <- c(body, word_section(tab$landscape, TRUE))
  }
  body <- c(body,word_section(tables[[length(tables)]]$landscape))
  write_xml('word/document.xml',paste0('<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>',paste0(body,collapse=''),'</w:body></w:document>'))
  write_xml('[Content_Types].xml','<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/><Override PartName="/word/footer.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/><Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/></Types>')
  write_xml('_rels/.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/></Relationships>')
  write_xml('word/_rels/document.xml.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/><Relationship Id="rIdFooter" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" Target="footer.xml"/></Relationships>')
  write_xml('word/styles.xml','<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/></w:rPr></w:rPrDefault></w:docDefaults><w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style><w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:after="120"/></w:pPr><w:rPr><w:b/><w:color w:val="000000"/><w:sz w:val="28"/></w:rPr></w:style><w:style w:type="paragraph" w:styleId="Note"><w:name w:val="Note"/><w:basedOn w:val="Normal"/><w:pPr><w:spacing w:before="120"/></w:pPr><w:rPr><w:sz w:val="18"/></w:rPr></w:style></w:styles>')
  write_xml('word/footer.xml','<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:fldSimple w:instr="PAGE"><w:r><w:rPr><w:sz w:val="18"/></w:rPr><w:t>1</w:t></w:r></w:fldSimple></w:p></w:ftr>')
  write_xml('docProps/core.xml','<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:title>AMPATH pregnancy and tuberculosis results</dc:title></cp:coreProperties>')
  files <- list.files(scratch,recursive=TRUE,all.files=TRUE,include.dirs=FALSE,no..=TRUE)
  Sys.setFileTime(file.path(scratch,files),as.POSIXct('2000-01-01',tz='UTC'))
  target <- normalizePath(output,mustWork=FALSE); old <- getwd(); on.exit(setwd(old),add=TRUE); setwd(scratch)
  if (file.exists(target)) unlink(target)
  status <- utils::zip(target,files=sort(files),flags='-q -X')
  if (status != 0L || !file.exists(target)) stop('Word archive creation failed.')
  invisible(target)
}
write_manuscript_html <- function(tables, output) {
  head <- '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AMPATH manuscript tables</title><style>body{font-family:Georgia,serif;margin:40px auto;max-width:1250px;padding:0 24px;color:#171717}h1{font-size:28px}h2{font-size:21px;margin-top:42px}table{border-collapse:collapse;width:100%;font-size:14px;border-top:2px solid;border-bottom:2px solid}th,td{padding:9px 10px;text-align:center;vertical-align:top}th{background:#f2f2f2;border-bottom:1px solid}td:first-child,th:first-child{text-align:left}.note{font-size:13px;line-height:1.5}section{overflow-x:auto;margin-bottom:40px}@media print{section{break-before:page}body{margin:0}}</style><body><h1>AMPATH pregnancy and tuberculosis manuscript tables</h1><p>Reconstructed cohort, SCCS, sensitivity and pathway results.</p>'
  parts <- vapply(tables,function(tab){
    rows <- vapply(seq_len(nrow(tab$data)),function(i)paste0('<tr>',paste0('<td>',xml_escape(as.character(tab$data[i,])),'</td>',collapse=''),'</tr>'),character(1))
    paste0('<section id="',tab$id,'"><h2>',gsub('\n','<br>',xml_escape(tab$title),fixed=TRUE),'</h2><table><thead><tr>',paste0('<th>',xml_escape(names(tab$data)),'</th>',collapse=''),'</tr></thead><tbody>',paste0(rows,collapse=''),'</tbody></table><p class="note">',xml_escape(tab$notes),'</p></section>')
  },character(1))
  writeLines(paste0(head,paste0(parts,collapse=''),'</body></html>'),output,useBytes=TRUE)
}
