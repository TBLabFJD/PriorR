# Administration of users and groups (only for administrators).

admin_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(7,
           shinydashboardPlus::box(
             title = "Users", width = NULL, status = "success", collapsible = FALSE,
             DT::DTOutput(ns("users")),
             hr(),
             uiOutput(ns("edit_user"))),
           shinydashboardPlus::box(
             title = "New user", width = NULL, status = "success", collapsible = TRUE, collapsed = TRUE,
             fluidRow(
               column(4, textInput(ns("new_username"), "User")),
               column(4, textInput(ns("new_full_name"), "Full name")),
               column(4, passwordInput(ns("new_password"), "Password (8+ characters)"))),
             fluidRow(
               column(4, checkboxInput(ns("new_is_admin"), "Administrator")),
               column(8, selectizeInput(ns("new_groups"), "Groups", choices = NULL, multiple = TRUE))),
             actionButton(ns("create_user"), "Create user", class = "btn-success"))),
    column(5,
           shinydashboardPlus::box(
             title = "Groups", width = NULL, status = "warning", collapsible = FALSE,
             helpText("Members of a group see the comments that the others share with their groups."),
             DT::DTOutput(ns("groups")),
             hr(),
             fluidRow(
               column(5, textInput(ns("new_group"), "New group")),
               column(7, textInput(ns("new_group_description"), "Description"))),
             actionButton(ns("create_group"), "Create group", class = "btn-warning"))))
}

admin_server <- function(id, user) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    changed <- reactiveVal(0)
    refresh <- function() changed(changed() + 1)
    is_admin <- reactive(isTRUE(user()$is_admin))

    users <- reactive({ req(is_admin()); changed(); list_users() })
    groups <- reactive({ req(is_admin()); changed(); list_groups() })

    observe(updateSelectizeInput(session, "new_groups", choices = setNames(groups()$id, groups()$name)))

    output$users <- DT::renderDT({
      u <- users()
      u$is_admin <- ifelse(u$is_admin, "yes", "")
      u$active <- ifelse(u$active, "yes", "no")
      DT::datatable(u[, c("username", "full_name", "groups", "is_admin", "active", "created_at")],
                    colnames = c("User", "Name", "Groups", "Administrator", "Active", "Created"),
                    rownames = FALSE, selection = "single", options = list(pageLength = 15, dom = "ftp"))
    })

    output$groups <- DT::renderDT({
      DT::datatable(groups()[, c("name", "description", "members")], colnames = c("Group", "Description", "Members"),
                    rownames = FALSE, selection = "none",
                    options = list(pageLength = 15, dom = "tp"))
    })

    selected_user <- reactive({
      i <- input$users_rows_selected
      req(length(i) == 1)
      users()[i, ]
    })

    output$edit_user <- renderUI({
      if (length(input$users_rows_selected) != 1) return(helpText("Select a user in the table to edit it."))
      u <- selected_user()
      member_of <- user_groups(u$id)$id
      tagList(
        h4(u$username),
        fluidRow(
          column(4, textInput(ns("edit_full_name"), "Full name", u$full_name)),
          column(4, checkboxInput(ns("edit_is_admin"), "Administrator", u$is_admin)),
          column(4, checkboxInput(ns("edit_active"), "Active (can sign in)", u$active))),
        selectizeInput(ns("edit_groups"), "Groups", choices = setNames(groups()$id, groups()$name),
                       selected = member_of, multiple = TRUE),
        fluidRow(
          column(6, passwordInput(ns("edit_password"), "New password (leave empty to keep it)")),
          column(6, br(), actionButton(ns("save_user"), paste("Save", u$username), class = "btn-success"))))
    })

    observeEvent(input$save_user, {
      req(is_admin())
      u <- selected_user()
      if (u$id == user()$id && (!isTRUE(input$edit_is_admin) || !isTRUE(input$edit_active)))
        return(showNotification("You can't remove your own administrator rights or deactivate yourself", type = "error"))
      result <- tryCatch({
        update_user(u$id, full_name = input$edit_full_name, is_admin = input$edit_is_admin, active = input$edit_active)
        set_user_groups(u$id, input$edit_groups)
        if (nzchar(input$edit_password)) set_password(u$id, input$edit_password)
      }, error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      showNotification(paste("Saved", u$username), type = "message")
      refresh()
    })

    observeEvent(input$create_user, {
      req(is_admin())
      result <- tryCatch({
        id <- create_user(input$new_username, input$new_password, input$new_full_name, input$new_is_admin)
        set_user_groups(id, input$new_groups)
      }, error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      showNotification(paste("User", input$new_username, "created"), type = "message")
      for (i in c("new_username", "new_full_name", "new_password")) updateTextInput(session, i, value = "")
      refresh()
    })

    observeEvent(input$create_group, {
      req(is_admin())
      result <- tryCatch(create_group(input$new_group, input$new_group_description), error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      updateTextInput(session, "new_group", value = "")
      updateTextInput(session, "new_group_description", value = "")
      refresh()
    })
  })
}
