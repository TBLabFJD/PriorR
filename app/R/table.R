# Variant tables: processed in the server (the browser only receives the visible page),
# fixed column widths, long values cut with "..." and shown in full on hover / click.

# Names shown for the columns: SYMBOL as GENE, spaces instead of "_", and sample columns
# (P001_GT...) starting with the field so they can be told apart when cut. The data keeps the
# original names (filters, profiles and downloads use them).
column_labels <- function(columns) {
  labels <- sub("^(.+)_(GT|VAF|AD|DP|SF|GD|VD|ROH)$", "\\2 \\1", columns)
  labels[columns == "SYMBOL"] <- "GENE"
  labels[columns == "GLOW_rank"] <- "GLOW rank"
  gsub("_", " ", labels, fixed = TRUE)
}

# Choices of a column picker: the labels shown, the original names as values
column_choices <- function(columns) stats::setNames(columns, column_labels(columns))

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

table_css <- shiny::tags$style(shiny::HTML("
  table.priorr-table { table-layout: fixed; }
  table.priorr-table th, table.priorr-table td { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  span.cell-full { cursor: zoom-in; }
"))

# DT table of variants. `row_color` (optional) colors whole rows; `full_cell_input` is the
# input that receives the clicked cell.
priorr_datatable <- function(DF, full_cell_input, row_color = NULL, page_length = 10) {
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
  column_defs <- c(
    list(list(targets = "_all", render = cell_render)),
    lapply(unique(widths), function(w) list(targets = which(widths == w) - 1L, width = w)),
    list(list(targets = which(widths == "70px") - 1L, className = "dt-center")),
    if (length(hidden)) list(list(targets = hidden, visible = FALSE, searchable = FALSE)))

  dt <- DT::datatable(
    DF, container = container, callback = cell_click_callback(full_cell_input),
    filter = "top", class = "display nowrap compact priorr-table", rownames = FALSE,
    selection = "single", escape = FALSE,
    options = list(scrollX = TRUE, autoWidth = FALSE, pageLength = page_length,
                   columnDefs = column_defs, searchHighlight = FALSE))
  if (!is.null(row_color)) {
    # light tones so that the text stays readable
    dt <- DT::formatStyle(dt, ".row_color", target = "row", backgroundColor = DT::styleEqual(
      c("salmon", "palegreen", "red", "darkorange"), c("#fbd3c4", "#c8f0c8", "#f7a1a1", "#fdcf94")))
  }
  dt
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
