# Comments and classifications of a variant, shared with the user's groups.
# Returns list(changed = reactive counter, open = function(variant) showing the dialog).
# A variant is a list(key, gene, label).

comments_ui <- function(id) {
  ns <- NS(id)
  tagList(uiOutput(ns("list")),
          actionButton(ns("add"), "Add comment / classification", icon = icon("comment"), class = "btn-sm"))
}

comments_server <- function(id, user, selected_variant) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    changed <- reactiveVal(0)
    dialog_variant <- reactiveVal(NULL)

    comments_of <- function(variant) {
      changed()
      visible_comments(user()$id, variant$key)
    }

    comments_html <- function(cm) {
      if (nrow(cm) == 0) return(p(class = "text-muted", "No comments visible for you."))
      tags$ul(class = "list-unstyled", lapply(seq_len(nrow(cm)), function(i) {
        c <- cm[i, ]
        shared <- switch(c$visibility, private = "only the author", groups = "author's groups",
                         group = paste("group", c$group_name))
        tags$li(style = "border-left: 3px solid #00a65a; padding-left: 8px; margin-bottom: 8px;",
                strong(c$author), span(class = "text-muted", sprintf(" %s · shared with %s", c$created_at, shared)),
                if (tolower(c$author) == tolower(user()$id)) tags$a(href = "#", icon("trash"), title = "Delete",
                                                    onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'}); return false;", ns("delete"), c$id)),
                if (!is.na(c$classification)) div(span(class = "label label-primary", c$classification),
                                                  if (nzchar(c$acmg %||% "")) span(class = "text-muted", style = "font-size: 12px;", " ACMG:", c$acmg)),
                if (!is.na(c$comment) && nzchar(c$comment)) div(style = "white-space: pre-wrap;", c$comment))
      }))
    }

    output$list <- renderUI({
      v <- selected_variant()
      req(v, user())
      comments_html(comments_of(v))
    })

    output$dialog_list <- renderUI({
      v <- dialog_variant()
      req(v)
      comments_html(comments_of(v))
    })

    open <- function(variant) {
      dialog_variant(variant)
      groups <- user_groups(user()$id)
      share <- c("All my groups" = "groups", setNames(sprintf("group:%s", groups$name), sprintf("Group %s", groups$name)),
                 "Only me" = "private")
      showModal(modalDialog(
        title = paste("Comments:", variant$label), size = "m", easyClose = TRUE,
        uiOutput(ns("dialog_list")),
        hr(),
        selectInput(ns("classification"), "Classification", c("(none)" = "", classifications)),
        selectizeInput(ns("acmg"), "ACMG criteria (suggested from the annotation: confirm, remove or add)",
                       choices = unique(c(variant$acmg, acmg_all_criteria)), selected = variant$acmg, multiple = TRUE, width = "100%"),
        uiOutput(ns("acmg_class")),
        textAreaInput(ns("comment"), "Comment", width = "100%", height = "90px"),
        selectInput(ns("share"), "Share with", share),
        footer = tagList(modalButton("Close"), actionButton(ns("save"), "Save", class = "btn-success"))))
    }

    output$acmg_class <- renderUI({
      cr <- input$acmg %||% character(0)
      helpText(sprintf("ACMG points %+d: %s", acmg_total(cr), acmg_class(cr)))
    })

    observeEvent(input$add, {
      v <- selected_variant()
      if (is.null(v)) return(showNotification("Select a variant in the table first", type = "warning"))
      open(v)
    })

    observeEvent(input$save, {
      v <- dialog_variant()
      req(v, user())
      share <- input$share
      group_id <- if (startsWith(share, "group:")) sub("^group:", "", share) else NA
      result <- tryCatch(
        add_comment(user()$id, v$key, gene = v$gene, classification = input$classification,
                    comment = input$comment, visibility = if (is.na(group_id)) share else "group", group_id = group_id,
                    acmg = input$acmg, snapshot = v$snapshot),
        error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      updateTextAreaInput(session, "comment", value = "")
      updateSelectInput(session, "classification", selected = "")
      changed(changed() + 1)
    })

    observeEvent(input$delete, {
      if (delete_comment(input$delete, user()$id)) changed(changed() + 1)
    })

    list(changed = changed, open = open)
  })
}
