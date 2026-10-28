# Analysis rounds of the SNV tab (rounds.R): per user and file. Runs inside snv_server;
# `params` is read when a round is saved. Returns the reactive with the variant keys to hide.

snv_rounds_server <- function(input, output, session, user, data, params, loaded) {
  ns <- session$ns
  loaded_file <- reactiveVal(NULL)
  rounds_changed <- reactiveVal(0)
  observeEvent(loaded(), loaded_file(loaded()$file))
  rounds <- reactive({
    req(user(), loaded_file())
    rounds_changed()
    list_rounds(user()$id, loaded_file())
  })
  hidden_keys <- reactive(if (isTRUE(input$hide_rounds) && !is.null(loaded_file())) unique(unlist(lapply(rounds(), `[[`, "keys"))))
  # Default name of a round: sample and main filters (panels, frequency, CADD...)
  round_default_name <- reactive({
    req(data())
    sample <- if (identical(input$family_source, "trio") && nzchar(input$trio_proband %||% "")) input$trio_proband
              else sample_names(names(data()))[1]
    catalog <- panel_catalog()
    panels <- catalog$name[match(c(input$panels, input$specific_panels), catalog$id)]
    af <- unlist(list(input$max_gnomadg_af, input$max_gnomadg_popmax, input$max_gnomade_af, input$max_fjd_af))
    af <- af[!is.na(af) & af < 1]
    round_name(sample = sample, panels = panels, hpo = length(input$phenotype), genes = split_genes(input$genes),
               max_af = if (length(af)) min(af), min_cadd = as.numeric(input$min_cadd %||% 0))
  })
  # A round only takes the genes (panels, HPO, genes), the frequencies and CADD: the other filters
  # (consequence, region, inheritance, scores...) are ignored, so the name says what it hides
  round_variants <- function() {
    p <- params()
    round_params <- snv_params(genes = p$genes, panel_genes = p$panel_genes, max_af = Filter(function(x) !is.null(x) && x < 1, p$max_af),
                               min_cadd = if (isTRUE(p$min_cadd > 0)) p$min_cadd)
    if (!length(round_params$genes) && is.null(round_params$panel_genes) && !length(round_params$max_af) && is.null(round_params$min_cadd))
      stop("A round needs gene panels / HPO / genes, a maximum frequency or a minimum CADD")
    filter_snv(data(), round_params)
  }
  observe(updateTextInput(session, "round_name", placeholder = round_default_name()))
  observeEvent(input$save_round, {
    req(user())
    if (is.null(loaded_file())) return(showNotification("Open a variant file first", type = "warning"))
    name <- if (nzchar(trimws(input$round_name %||% ""))) input$round_name else round_default_name()
    DF <- NULL
    result <- tryCatch({
      DF <- round_variants()
      add_round(user()$id, loaded_file(), variant_key(DF$CHROM, DF$POS, DF$REF, DF$ALT), name)
    }, error = function(e) e)
    if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
    updateTextInput(session, "round_name", value = "")
    rounds_changed(rounds_changed() + 1)
    showNotification(sprintf("Round saved: %d variants. Tick 'Hide the variants of previous rounds' and change the filters.",
                             length(unique(variant_key(DF$CHROM, DF$POS, DF$REF, DF$ALT)))), type = "message")
  })
  observeEvent(input$delete_round, {
    delete_round(user()$id, loaded_file(), as.integer(input$delete_round))
    rounds_changed(rounds_changed() + 1)
  })
  output$rounds <- renderUI({
    r <- rounds()
    if (!length(r)) return(div(style = "margin: 6px 15px; font-size: 12px; white-space: normal;", "No rounds saved for this file."))
    div(style = "margin: 6px 15px; font-size: 12px; white-space: normal;",
        lapply(seq_along(r), function(i) div(style = "margin-bottom: 4px;",
          icon("layer-group"), strong(r[[i]]$name), sprintf(" · %d variants · %s", length(r[[i]]$keys), substr(r[[i]]$created_at, 1, 16)),
          if (nzchar(r[[i]]$summary %||% "")) span(style = "color: #b8c7ce;", paste0(" (", r[[i]]$summary, ")")),
          tags$a(href = "#", style = "margin-left: 6px; color: #ef9a9a;", title = "Delete the round", icon("trash"),
                 onclick = sprintf("Shiny.setInputValue('%s', %d, {priority: 'event'}); return false;", ns("delete_round"), i)))))
  })
  hidden_keys
}
