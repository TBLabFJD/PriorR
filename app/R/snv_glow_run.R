# Run GLOWgenes with Docker in the background from the SNV tab. The run works in the shared folder
# (glowgenes/work/<name>), so it survives a restart of PriorR and the whole lab sees it;
# the ranking goes to glowgenes/<name>.tsv (glowgenes.R).

snv_glow_run_server <- function(input, output, session, user) {
  ns <- session$ns
  output$glow_docker <- renderUI({
    problem <- docker_problem()
    if (nzchar(problem))
      div(style = "margin: 4px 15px; font-size: 12px; color: #ffcc80; white-space: normal;", icon("exclamation-triangle"), problem,
          br(), "Meanwhile, use the precomputed GLOWgenes rankings.")
  })
  my_runs <- reactiveVal(character(0))
  glow_runs <- reactivePoll(5000, session,
    checkFunc = function() file.mtime(Sys.glob(file.path(glow_runs_dir(), "work", "*", "status.json"))),
    valueFunc = glow_run_statuses)
  observeEvent(input$glow_run, {
    name <- gsub("[^A-Za-z0-9 _.-]", "", trimws(input$glow_run_name))
    genes <- if (!is.null(input$glow_run_file)) read_gene_list(input$glow_run_file$datapath) else panel_genes(input$glow_run_panel)
    if (!nzchar(name)) return(showNotification("Give the ranking a name", type = "error"))
    if (length(genes) < glow_min_genes) return(showNotification(sprintf("Choose a panel or a file with at least %d genes", glow_min_genes), type = "error"))
    if (file.exists(file.path(glow_runs_dir(), paste0(name, ".tsv")))) return(showNotification("There is already a ranking with that name", type = "error"))
    if (identical(glow_runs()$stage[glow_runs()$name == name], "running")) return(showNotification("That ranking is already running", type = "error"))
    problem <- docker_problem()
    if (nzchar(problem)) return(showNotification(problem, type = "error", duration = NULL))
    workdir <- file.path(glow_runs_dir(), "work", name)
    unlink(workdir, recursive = TRUE)   # a failed run with the same name
    dir.create(workdir, recursive = TRUE)
    writeLines(genes, file.path(workdir, "genes.txt"))
    jsonlite::write_json(list(stage = "running", message = "Starting Docker...", user = user()$username, time = format(Sys.time())),
                         file.path(workdir, "status.json"), auto_unbox = TRUE)
    system2(file.path(R.home("bin"), "Rscript"),
            c(shQuote(normalizePath(file.path(data_dir, "scripts", "run_glowgenes.R"))), shQuote(normalizePath(workdir)),
              shQuote(file.path(normalizePath(glow_runs_dir()), paste0(name, ".tsv"))), glow_docker_image), wait = FALSE)
    my_runs(c(my_runs(), name))
  })
  # When a run of this session finishes, its ranking is selected
  finished <- reactive({
    req(length(my_runs()) > 0)
    glow_runs()
    my_runs()[file.exists(file.path(glow_runs_dir(), paste0(my_runs(), ".tsv")))]
  })
  observeEvent(finished(), {
    req(length(finished()) > 0)
    updateSelectizeInput(session, "glow_panels", choices = glow_rank_choices(),
                         selected = unique(c(isolate(input$glow_panels), paste0("run:", finished()))))
    updateRadioButtons(session, "glow_source", selected = "panel")
    showNotification(paste("GLOWgenes ranking ready:", paste(finished(), collapse = ", ")), type = "message")
    my_runs(setdiff(my_runs(), finished()))
  })
  output$glow_run_status <- renderUI({
    runs <- glow_runs()
    req(nrow(runs) > 0)
    div(style = "margin: 6px 15px; font-size: 12px; white-space: normal;",
        lapply(seq_len(nrow(runs)), function(i) div(
          if (runs$stage[i] == "running") tagList(icon("spinner", class = "fa-spin"), strong(runs$name[i]), sprintf("(%s, since %s): %s", runs$user[i], runs$time[i], runs$message[i]))
          else span(style = "color: #ef9a9a;", icon("times"), strong(runs$name[i]), runs$message[i]))))
  })
}
