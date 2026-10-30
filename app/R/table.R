# Variant tables: processed in the server (the browser only receives the visible page),
# fixed column widths, long values cut with "..." and shown in full on hover / click.

# Names shown for the columns: SYMBOL as GENE, spaces instead of "_", and sample columns
# (P001_GT...) starting with the field so they can be told apart when cut. The data keeps the
# original names (filters, profiles and downloads use them).
column_labels <- function(columns) {
  labels <- sub("^(.+)_(GT|VAF|AD|DP|SF|GD|VD|ROH)$", "\\2 \\1", columns)
  labels[columns == "SYMBOL"] <- "GENE"
  labels[columns == "GLOW_rank"] <- "GLOW rank"
  labels[columns == "Variants_in_gene"] <- "Var. in gene"
  gsub("_", " ", labels, fixed = TRUE)
}

# Choices of a column picker: the labels shown, the original names as values
column_choices <- function(columns) stats::setNames(columns, column_labels(columns))

# Choices of the column picker: the displayed columns first, in their order (the order of the
# picker is the order of the table), then the others
ordered_column_choices <- function(columns, shown) column_choices(unique(c(intersect(shown, columns), columns)))

# Columns in the order of the headers after dragging one (ColReorder): `order` (column names) moves
# the shown columns, the selected ones that are not shown stay at the end
reordered_columns <- function(selected, order) c(intersect(order, selected), setdiff(selected, order))

column_width <- function(columns) {
  narrow <- grepl("^(CHROM|REF|ALT|CANONICAL|STRAND|EXON|INTRON|Comments|GLOW_rank)$|_(GT|VAF|AD|DP|SF|GD|VD)$", columns)
  wide <- columns %in% c("HGVSc", "HGVSp", "Consequence", "CLNSIG", "Gene_full_name", "OMIM_phenotype",
                         "Orphanet_disorder", "Existing_variation", "Classification")
  ifelse(narrow, "70px", ifelse(wide, "170px", "110px"))
}

# Renders long values in a span with the full value as tooltip (cut at 200 characters in the page)
cell_render <- DT::JS("function(data, type, row, meta) {
  if (type !== 'display' || data === null || data === undefined) return data;
  var s = String(data);
  if (s.length <= 8 || s.charAt(0) === '<') return s;
  var esc = function(t) { return t.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/\"/g, '&quot;'); };
  var tip = s.length > 1000 ? s.substr(0, 1000) + '\\u2026 (click to see all)' : s;
  var shown = s.length > 200 ? s.substr(0, 200) + '\\u2026' : s;
  return '<span class=\"cell-full\" title=\"' + esc(tip) + '\">' + esc(shown) + '</span>';
}")

# Click on a cut value -> input `full_cell_input` = {column, text}
cell_click_callback <- function(full_cell_input) {
  DT::JS(sprintf("table.on('click', 'span.cell-full', function() {
    var td = $(this).closest('td')[0];
    var cell = table.cell(td);
    if (td.scrollWidth <= td.clientWidth && String(cell.data()).length <= 200) return;
    var header = table.column(cell.index().column).header();
    Shiny.setInputValue('%s', {column: $(header).attr('title') || $(header).text(), text: String(cell.data())}, {priority: 'event'});
  });", full_cell_input))
}

# After dragging a header: the column names (title of the headers) in the new order
column_reorder_callback <- function(input_id) {
  sprintf("table.on('column-reorder', function() {
    var cols = [];
    table.columns().every(function() { var h = $(this.header()).attr('title'); if (h) cols.push(h); });
    Shiny.setInputValue('%s', cols, {priority: 'event'});
  });", input_id)
}

# Small colored labels of some columns (ClinVar, genotypes, ACMG class of CNVs, gene in bold)
esc_js <- "var esc = function(t) { return String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/\"/g, '&quot;'); };"
render_gene <- DT::JS(paste("function(data, type) {", esc_js,
  "if (type !== 'display' || data === null || data === undefined || data === '') return data;",
  "return '<span class=\"pr-gene-cell\">' + esc(data) + '</span>'; }"))
render_clinvar <- DT::JS(paste("function(data, type) {", esc_js,
  "if (type !== 'display' || data === null || data === undefined || data === '') return data;",
  "var s = String(data).toLowerCase().replace(/ /g, '_'), c = 'other';",
  "if (s.indexOf('conflicting') >= 0) c = 'conflicting';",
  "else if (s.indexOf('likely_pathogenic') >= 0 && s.indexOf('pathogenic') !== 0) c = 'lp';",
  "else if (s.indexOf('pathogenic') >= 0) c = 'p';",
  "else if (s.indexOf('uncertain') >= 0) c = 'vus';",
  "else if (s.indexOf('benign') >= 0) c = 'benign';",
  "return '<span class=\"pr-pill pr-cv-' + c + '\" title=\"' + esc(data) + '\">' + esc(String(data).replace(/_/g, ' ')) + '</span>'; }"))
render_genotype <- DT::JS(paste("function(data, type) {", esc_js,
  "if (type !== 'display' || data === null || data === undefined || data === '') return data;",
  "var g = String(data), c = /^(1[\\/|]1|1)$/.test(g) ? 'hom' : /^(0[\\/|]1|1[\\/|]0)$/.test(g) ? 'het' : /^(0[\\/|]0|0)$/.test(g) ? 'ref' : 'other';",
  "return '<span class=\"pr-gt pr-gt-' + c + '\">' + esc(g) + '</span>'; }"))
render_acmg_class <- DT::JS(paste("function(data, type) {",
  "if (type !== 'display' || data === null || data === undefined || data === '') return data;",
  "var n = parseInt(data, 10); if (isNaN(n)) return data;",
  "return '<span class=\"pr-pill pr-acmg-' + n + '\">' + n + '</span>'; }"))

table_css <- shiny::tags$style(shiny::HTML("
  table.priorr-table { table-layout: fixed; }
  table.priorr-table th, table.priorr-table td { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  span.cell-full { cursor: zoom-in; }
  /* header and rows */
  table.priorr-table thead th { background: #f5f7f8; color: #37474f; border-bottom: 2px solid #00a65a !important; }
  table.priorr-table tbody tr:hover td { background-color: #eef6f1 !important; }
  table.priorr-table tbody td { border-top: 1px solid #eef1f3; }
  .pr-gene-cell { font-weight: 700; color: #1b5e20; }
  /* colored rows: a light tint and a bar on the left (priority, ROH, not in the other test) */
  table.priorr-table tbody tr.pr-row-red td { background-color: #fdeceb; }
  table.priorr-table tbody tr.pr-row-red td:first-child { box-shadow: inset 5px 0 0 #d32f2f; }
  table.priorr-table tbody tr.pr-row-orange td { background-color: #fff3e3; }
  table.priorr-table tbody tr.pr-row-orange td:first-child { box-shadow: inset 5px 0 0 #f57c00; }
  table.priorr-table tbody tr.pr-row-roh td { background-color: #fbf1ec; }
  table.priorr-table tbody tr.pr-row-roh td:first-child { box-shadow: inset 5px 0 0 #a1887f; }
  table.priorr-table tbody tr.pr-row-new td { background-color: #edf7ee; }
  table.priorr-table tbody tr.pr-row-new td:first-child { box-shadow: inset 5px 0 0 #43a047; }
  /* selected row: a soft blue instead of the strong one of DataTables; the colored bar stays */
  table.dataTable.display.priorr-table > tbody > tr.selected > td, table.dataTable.display.priorr-table > tbody > tr.selected:hover > td {
    background-color: #d6e8fb !important; box-shadow: none !important; color: #0d2a4a !important; }
  table.dataTable.display.priorr-table > tbody > tr.selected > td .pr-gene-cell { color: #0d47a1; }
  table.dataTable.display.priorr-table > tbody > tr.selected.pr-row-red > td:first-child { box-shadow: inset 5px 0 0 #d32f2f !important; }
  table.dataTable.display.priorr-table > tbody > tr.selected.pr-row-orange > td:first-child { box-shadow: inset 5px 0 0 #f57c00 !important; }
  table.dataTable.display.priorr-table > tbody > tr.selected.pr-row-roh > td:first-child { box-shadow: inset 5px 0 0 #a1887f !important; }
  table.dataTable.display.priorr-table > tbody > tr.selected.pr-row-new > td:first-child { box-shadow: inset 5px 0 0 #43a047 !important; }
  /* labels inside cells */
  .pr-pill { display: inline-block; padding: 0 7px; border-radius: 9px; font-size: 11.5px; font-weight: 600; line-height: 17px;
             max-width: 100%; overflow: hidden; text-overflow: ellipsis; vertical-align: middle; }
  .pr-cv-p { background: #ffcdd2; color: #b71c1c; } .pr-cv-lp { background: #ffe0b2; color: #e65100; }
  .pr-cv-vus { background: #fff8c5; color: #8d6e00; } .pr-cv-benign { background: #dcedc8; color: #33691e; }
  .pr-cv-conflicting { background: #e1bee7; color: #6a1b9a; } .pr-cv-other { background: #eceff1; color: #546e7a; }
  .pr-acmg-1 { background: #c8e6c9; color: #1b5e20; } .pr-acmg-2 { background: #e6f4c8; color: #558b2f; }
  .pr-acmg-3 { background: #fff3b0; color: #8d6e00; } .pr-acmg-4 { background: #fdcf94; color: #bf360c; }
  .pr-acmg-5 { background: #f7a1a1; color: #b71c1c; }
  .pr-gt { display: inline-block; min-width: 30px; padding: 0 5px; border-radius: 8px; font-size: 11.5px; font-weight: 600; line-height: 16px; }
  .pr-gt-het { background: #e3f2fd; color: #1565c0; } .pr-gt-hom { background: #ede7f6; color: #5e35b1; }
  .pr-gt-ref { background: #f1f8e9; color: #558b2f; } .pr-gt-other { color: #90a4ae; }
"))

# CSS class of each row color (row_colors)
row_classes <- c(red = "pr-row-red", darkorange = "pr-row-orange", salmon = "pr-row-roh", palegreen = "pr-row-new")

# DT table of variants. `row_color` (optional) colors whole rows; `full_cell_input` is the
# input that receives the clicked cell.
# Allele frequency columns (gnomAD, FJD, CSVS...; not the VAF of the samples): empty = not seen
# in the population, so it sorts as 0
is_af_column <- function(columns) grepl("(^|_)(AF|MAF)(_|$)", columns) & !grepl("_VAF$", columns)

# Sort key of a column: empty frequencies as 0; text columns with numbers sort as numbers
sort_key <- function(x, name) {
  if (!is.numeric(x)) {
    num <- suppressWarnings(as.numeric(as.character(x)))
    if (all(is.na(num) == is.na(x))) x <- num
  }
  if (is.numeric(x)) {
    if (is_af_column(name)) x[is.na(x)] <- 0
    return(x)
  }
  xtfrm(x)
}

# Server side processing of the tables (DT::renderDT funcFilter): DT's search and pages, with
# the order of sort_key. The rows of DT_rows_all / DT_rows_current stay those of `data`, so the
# selected rows are right.
priorr_dt_filter <- function(data, params) {
  dt_filter <- utils::getFromNamespace("dataTablesFilter", "DT")
  keys <- list()
  for (o in params$order) {
    col <- params$columns[[o[["column"]]]]
    if (!identical(col[["orderable"]], "true")) next
    j <- match(col[["name"]], names(data))
    if (is.na(j)) j <- as.integer(o[["column"]]) + 1L
    key <- sort_key(data[[j]], names(data)[j])
    keys[[length(keys) + 1]] <- if (identical(o[["dir"]], "desc")) -key else key
  }
  if (!length(keys)) return(dt_filter(data, params))
  perm <- do.call(order, c(keys, list(na.last = TRUE)))
  params$order <- list()
  res <- dt_filter(data[perm, , drop = FALSE], params)
  res$DT_rows_all <- perm[res$DT_rows_all]
  res$DT_rows_current <- perm[res$DT_rows_current]
  res
}

# Sorting by a column keeps the previous order within equal values: the columns clicked before
# go after it in the order sent to the server (shift + click still chooses several columns)
sort_history_callback <- DT::JS("var sortHistory = [];
  table.on('preXhr.dt', function(e, settings, d) {
    if (!d.order || !d.order.length) return;
    if (d.order.length === 1) {
      var o = d.order[0];
      sortHistory = [o].concat(sortHistory.filter(function(h) { return h.column != o.column; })).slice(0, 5);
      d.order = sortHistory.slice();
    } else sortHistory = d.order.slice();
  });")

# reorder_input: dragging a header sets this input to the column names in their new order
priorr_datatable <- function(DF, full_cell_input, row_color = NULL, page_length = 10, reorder_input = NULL) {
  hidden <- integer(0)
  if (!is.null(row_color)) {
    DF$.row_color <- row_color
    hidden <- ncol(DF) - 1L
  }
  headers <- names(DF)
  labels <- column_labels(headers)
  labels[headers == ".row_color"] <- ""
  container <- tags$table(class = "display priorr-table", tags$thead(tags$tr(
    mapply(function(h, l) tags$th(title = h, l), headers, labels, SIMPLIFY = FALSE, USE.NAMES = FALSE))))
  widths <- column_width(headers)
  # specific renderers first: the first definition of a column wins
  target <- function(x) which(x) - 1L
  special <- Filter(function(d) length(d$targets) > 0, list(
    list(targets = target(headers == "SYMBOL"), render = render_gene),
    list(targets = target(headers %in% c("CLNSIG", "Classification")), render = render_clinvar),
    list(targets = target(grepl("_GT$", headers)), render = render_genotype),
    list(targets = target(headers == "ACMG_class"), render = render_acmg_class)))
  column_defs <- c(
    special,
    list(list(targets = "_all", render = cell_render)),
    lapply(unique(widths), function(w) list(targets = which(widths == w) - 1L, width = w)),
    list(list(targets = which(widths == "70px" | headers == "Variants_in_gene") - 1L, className = "dt-center")),
    if (length(hidden)) list(list(targets = hidden, visible = FALSE, searchable = FALSE)))

  # colored rows: a CSS class from the hidden .row_color column (table_css)
  row_callback <- if (length(hidden)) DT::JS(sprintf(
    "function(row, data) { var c = %s[data[%d]]; $(row).removeClass('pr-row-red pr-row-orange pr-row-roh pr-row-new'); if (c) $(row).addClass(c); }",
    jsonlite::toJSON(as.list(row_classes), auto_unbox = TRUE), hidden))
  options <- list(scrollX = TRUE, autoWidth = FALSE, pageLength = page_length,
                  columnDefs = column_defs, searchHighlight = FALSE)
  if (!is.null(row_callback)) options$rowCallback <- row_callback
  callback <- DT::JS(paste(cell_click_callback(full_cell_input), sort_history_callback, sep = "\n"))
  if (!is.null(reorder_input)) {
    options$colReorder <- list(realtime = FALSE)
    callback <- DT::JS(paste(callback, column_reorder_callback(reorder_input), sep = "\n"))
  }
  DT::datatable(
    DF, container = container, callback = callback,
    filter = "top", class = "display nowrap compact priorr-table", rownames = FALSE,
    selection = "single", escape = FALSE, options = options,
    extensions = if (!is.null(reorder_input)) "ColReorder" else list())
}

show_full_cell <- function(cell) {
  text <- gsub("(;|&)\\s*", "\\1\n", cell$text)
  showModal(modalDialog(
    title = cell$column, size = "l", easyClose = TRUE, footer = modalButton("Close"),
    tags$pre(style = "white-space: pre-wrap; word-break: break-word; max-height: 65vh; overflow-y: auto;", text)))
}

# Row colors: ROH (salmon) < not in the other test (green) < prioritization (orange, red)
row_colors <- function(DF, roh = FALSE, other_test = FALSE, priority = FALSE) {
  color <- rep("", nrow(DF))
  if (roh) color[attr(DF, "roh") %in% TRUE] <- "salmon"
  if (other_test) color[attr(DF, "not_in_other_test") %in% TRUE] <- "palegreen"
  if (priority) {
    color[attr(DF, "priority") %in% "orange"] <- "darkorange"
    color[attr(DF, "priority") %in% "red"] <- "red"
  }
  color
}

# A "Classify" button per row; clicking it sets `input_id` to the row number
classify_buttons <- function(n, input_id) {
  sprintf('<button class="btn btn-default btn-xs" onclick="Shiny.setInputValue(\'%s\', %d, {priority: \'event\'})">Classify</button>',
          input_id, seq_len(n))
}

# Bootstrap modal that stays in the page (its inputs keep their values when closed)
static_modal <- function(id, title, ..., size = "lg") {
  div(class = "modal fade", id = id, tabindex = "-1", role = "dialog",
      div(class = paste0("modal-dialog modal-", size), role = "document",
          div(class = "modal-content",
              div(class = "modal-header",
                  tags$button(type = "button", class = "close", `data-dismiss` = "modal", HTML("&times;")),
                  h4(class = "modal-title", title)),
              div(class = "modal-body", ...))))
}

toggle_modal <- function(id, show = TRUE) {
  shinyjs::runjs(sprintf("$('#%s').modal('%s');", id, if (show) "show" else "hide"))
}
