# Reads of the selected variant with igv.js (the web version of IGV). The BAM / CRAM files stay in
# the user's computer: PriorR serves only the pieces igv.js asks for (HTTP range requests), and only
# of the files chosen in the session. The reference genome and the gene track come from the IGV servers.

igv_js_url <- "https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"

reads_ui <- function(id) {
  ns <- NS(id)
  shinydashboardPlus::box(
    title = "Reads (IGV)", width = NULL, status = "primary", collapsible = TRUE,
    singleton(tags$head(tags$script(src = igv_js_url), tags$script(HTML(igv_handler_js)))),
    fluidRow(
      column(8,
             shinyFiles::shinyFilesButton(ns("files"), "Open BAM / CRAM files", "Open the BAM or CRAM files (with their .bai / .crai index)",
                                          multiple = TRUE, icon = icon("folder-open")),
             actionButton(ns("igv_desktop"), "Open in IGV desktop", icon = icon("external-link-alt")),
             actionLink(ns("clear"), "Close the BAM / CRAM files", style = "margin-left: 10px;",
                        title = "Takes the files out of the viewer (nothing is deleted from the disk)")),
      column(4, uiOutput(ns("files_info")))),
    uiOutput(ns("hint")),
    div(id = ns("igv"), style = "margin-top: 10px;"))
}

# Creates the IGV browser once and moves it to each variant; recreated when the files change
igv_handler_js <- "
Shiny.addCustomMessageHandler('priorr-igv', async function(m) {
  var el = document.getElementById(m.id);
  if (!el || typeof igv === 'undefined') return;
  var state = window.priorrIgv || {};
  if (state.browser && state.key === m.key) { state.browser.search(m.locus); return; }
  if (state.browser) { igv.removeBrowser(state.browser); }
  el.innerHTML = '';
  var browser = await igv.createBrowser(el, {genome: m.genome, locus: m.locus, tracks: m.tracks});
  window.priorrIgv = {browser: browser, key: m.key};
});
Shiny.addCustomMessageHandler('priorr-igv-desktop', async function(urls) {
  for (const u of urls) { try { await fetch(u, {mode: 'no-cors'}); } catch (e) {} }
});
"

is_alignment_file <- function(path) grepl("\\.(bam|cram)$", path, ignore.case = TRUE)

# Index of a BAM / CRAM file next to it (file.bam.bai, file.bai, file.bam.csi, file.cram.crai, file.crai;
# also in capitals)
alignment_index <- function(path) {
  ext <- if (grepl("\\.cram$", path, ignore.case = TRUE)) "crai" else c("bai", "csi")
  ext <- c(ext, toupper(ext))
  base <- tools::file_path_sans_ext(path)
  candidates <- c(paste0(path, ".", ext), paste0(base, ".", ext))
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1] else NA_character_
}

# Sample of a file: the longest sample name contained in the file name
file_sample <- function(path, samples) {
  hits <- samples[vapply(samples, function(s) grepl(s, basename(path), fixed = TRUE), logical(1))]
  if (length(hits)) hits[which.max(nchar(hits))] else tools::file_path_sans_ext(basename(path))
}

# Answers igv.js requests of a piece of a file ("Range: bytes=start-end")
range_response <- function(data, req) {
  size <- file.info(data$path)$size
  range <- req$HTTP_RANGE
  headers <- list(`Accept-Ranges` = "bytes", `Cache-Control` = "no-cache")
  if (is.null(range) || !nzchar(range)) {
    if (size > 200e6) return(shiny::httpResponse(416L, "text/plain", "Range requests only", headers))
    return(shiny::httpResponse(200L, "application/octet-stream", readBin(data$path, "raw", size), headers))
  }
  bounds <- as.numeric(strsplit(sub("^bytes=", "", range), "-", fixed = TRUE)[[1]])
  from <- bounds[1]
  to <- if (length(bounds) > 1 && !is.na(bounds[2])) min(bounds[2], size - 1) else size - 1
  if (is.na(from) || from >= size) {
    return(shiny::httpResponse(416L, "text/plain", "", c(headers, `Content-Range` = sprintf("bytes */%.0f", size))))
  }
  con <- file(data$path, "rb")
  on.exit(close(con))
  seek(con, from)
  content <- readBin(con, "raw", to - from + 1)
  shiny::httpResponse(206L, "application/octet-stream", content,
                      c(headers, `Content-Range` = sprintf("bytes %.0f-%.0f/%.0f", from, to, size)))
}

# variant: reactive with the selected row; data: the loaded variants (to know the samples)
# Region shown for a SNV / indel: the variant +- 60 bp
snv_locus <- function(v) {
  values <- link_values(v)
  if (is.na(values$pos)) return(NULL)
  sprintf("%s:%.0f-%.0f", values$chrom, as.numeric(values$pos) - 60, as.numeric(values$pos) + 60)
}

# Region shown for a CNV: the whole CNV when it is small, its two breakpoints when it is big
# (reads are only loaded for small windows)
cnv_locus <- function(v) {
  chrom <- as.character(v$SV_chrom %||% NA)
  start <- suppressWarnings(as.numeric(v$SV_start %||% NA)); end <- suppressWarnings(as.numeric(v$SV_end %||% NA))
  if (is.na(chrom) || is.na(start) || is.na(end)) return(NULL)
  if (end - start <= 20000) return(sprintf("%s:%.0f-%.0f", chrom, max(1, start - 500), end + 500))
  sprintf("%s:%.0f-%.0f %s:%.0f-%.0f", chrom, max(1, start - 300), start + 300, chrom, max(1, end - 300), end + 300)
}

reads_server <- function(id, variant, data, locus_of = snv_locus) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    roots <- file_roots()
    # only BAM / CRAM files are listed (also in capitals); their indexes are found automatically
    shinyFiles::shinyFileChoose(input, "files", roots = roots, filetypes = c("bam", "cram", "BAM", "CRAM"))
    files <- reactiveVal(data.frame(sample = character(0), path = character(0), index = character(0),
                                    format = character(0), stringsAsFactors = FALSE))

    add_files <- function(paths) {
      other <- paths[!is_alignment_file(paths)]
      if (length(other)) showNotification(paste("Not BAM / CRAM files (the indexes are found automatically):",
                                                paste(basename(other), collapse = ", ")), type = "warning")
      paths <- paths[is_alignment_file(paths)]
      if (!length(paths)) return()
      samples <- if (is.null(data())) character(0) else sample_names(names(data()))
      new <- data.frame(path = paths, stringsAsFactors = FALSE)
      new$index <- vapply(paths, alignment_index, "")
      new$sample <- vapply(paths, file_sample, "", samples = samples)
      new$format <- ifelse(grepl("\\.cram$", paths, ignore.case = TRUE), "cram", "bam")
      missing <- new$path[is.na(new$index)]
      if (length(missing)) showNotification(paste("Without index (.bai / .crai), not shown:", paste(basename(missing), collapse = ", ")),
                                            type = "error", duration = NULL)
      new <- new[!is.na(new$index), ]
      all <- rbind(files(), new[, names(files())])
      files(all[!duplicated(all$path), ])
    }

    observeEvent(input$files, {
      chosen <- shinyFiles::parseFilePaths(roots, input$files)
      if (nrow(chosen)) add_files(as.character(chosen$datapath))
    })
    observeEvent(input$clear, files(files()[0, ]))

    output$files_info <- renderUI({
      f <- files()
      if (!nrow(f)) return(NULL)
      tags$ul(style = "padding-left: 18px; margin: 0;", lapply(seq_len(nrow(f)), function(i)
        tags$li(strong(f$sample[i]), ": ", basename(f$path[i]))))
    })

    output$hint <- renderUI({
      if (!nrow(files())) helpText("Open the BAM or CRAM files of the samples (with their index) to see the reads of the selected variant.",
                                   "The files are read from this computer; the reference genome comes from the IGV servers (internet).")
      else if (is.null(variant())) helpText("Select a variant in the table.")
    })

    # Each file and its index get a private URL of this session
    tracks <- reactive({
      f <- files()
      lapply(seq_len(nrow(f)), function(i) {
        url <- session$registerDataObj(sprintf("aln%d", i), list(path = f$path[i]), range_response)
        index <- session$registerDataObj(sprintf("aln%d_index", i), list(path = f$index[i]), range_response)
        list(name = f$sample[i], type = "alignment", format = f$format[i], url = url, indexURL = index,
             height = 300, displayMode = "EXPANDED")
      })
    })

    locus <- reactive(req(locus_of(req(variant()))))

    observe({
      req(nrow(files()) > 0)
      session$sendCustomMessage("priorr-igv", list(
        id = ns("igv"), genome = if (genome_build() == "GRCh37") "hg19" else "hg38", locus = locus(), tracks = tracks(),
        key = paste(files()$path, collapse = "|")))
    })

    # IGV desktop (batch port 60151): load the files and go to the variant
    observeEvent(input$igv_desktop, {
      f <- files()
      if (!nrow(f) || is.null(variant())) return(showNotification("Open the BAM / CRAM files and select a variant first", type = "warning"))
      base <- "http://127.0.0.1:60151"
      urls <- c(sprintf("%s/load?file=%s&merge=true", base, vapply(f$path, utils::URLencode, "", reserved = TRUE)),
                sprintf("%s/goto?locus=%s", base, utils::URLencode(locus(), reserved = TRUE)))
      session$sendCustomMessage("priorr-igv-desktop", as.list(unname(urls)))
      showNotification("Sent to IGV desktop (IGV must be open, with 'Enable port' in View > Preferences > Advanced)", type = "message")
    })

    # Saved sessions remember the files
    onBookmark(function(state) state$values$alignment_files <- files()$path)
    onRestored(function(state) {
      paths <- unlist(state$values$alignment_files)
      if (length(paths)) isolate(add_files(paths[file.exists(paths)]))
    })
  })
}
