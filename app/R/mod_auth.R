# Login, logout and password change. Returns a reactive with the signed-in user
# (list with id, username, full_name, is_admin) or NULL.

auth_header_ui <- function(id) {
  ns <- NS(id)
  tags$li(class = "dropdown", style = "padding: 8px 15px;", uiOutput(ns("user_menu"), inline = TRUE))
}

auth_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    user <- reactiveVal(NULL)

    login_dialog <- function(message = NULL) {
      first_run <- user_count() == 0
      remember <- if (remember_enabled())
        checkboxInput(ns("remember"), "Remember me on this computer (don't ask for the password again)", value = TRUE)
      modalDialog(
        title = if (first_run) paste("Welcome to", priorr_name) else "Sign in",
        easyClose = FALSE, footer = NULL, size = "s",
        div(style = "text-align: center; margin-bottom: 15px;",
            tags$img(src = "logo_v4.png", style = "width: 220px;"),
            div(class = "text-muted", paste("Version", priorr_version))),
        if (first_run) div(class = "alert alert-info", style = "font-size: 13px;",
          p("This is the first time PriorR is opened with this shared folder."),
          p("Create your user: you will be the ", strong("administrator"), " and you will add your colleagues and ",
            "create groups in the 'Users and groups' tab. This is done only once.")),
        textInput(ns("username"), "User"),
        if (first_run) textInput(ns("full_name"), "Full name"),
        passwordInput(ns("password"), if (first_run) "Password (8 characters or more)" else "Password"),
        if (first_run) passwordInput(ns("password2"), "Repeat the password"),
        remember,
        if (!is.null(message)) div(class = "text-danger", style = "margin-bottom: 10px;", message),
        actionButton(ns(if (first_run) "create_admin" else "login"), if (first_run) "Create and sign in" else "Sign in",
                     class = "btn-success", style = "width: 100%;"),
        if (!first_run) div(class = "text-muted", style = "font-size: 11.5px; margin-top: 10px;",
          "No user yet? An administrator creates it in 'Users and groups'. Forgotten password of the administrator: see Help > Users and groups."),
        div(style = "text-align: center; margin-top: 15px; font-size: 11px; color: #90a4ae;",
            tags$a(href = tblab_url, target = "_blank", style = "color: #90a4ae;",
                   tags$img(src = "tblab_logo.png", style = "width: 26px; vertical-align: middle; margin-right: 4px;"),
                   "Bioinformatics Unit \u00b7 IIS-Fundación Jiménez Díaz")),
        tags$script(sprintf("$(document).on('keyup', '#%s', function(e) { if (e.key === 'Enter') $('#%s').click(); });",
                            ns("password"), ns("login"))))
    }

    signed_in <- function(u) {
      if (remember_enabled() && isTRUE(input$remember)) remember_user(u$id)
      user(u)
      removeModal()
    }

    # A computer that remembers its user goes straight in
    remembered <- remembered_user()
    if (is.null(remembered)) showModal(login_dialog()) else user(remembered)

    observeEvent(input$login, {
      u <- authenticate(input$username, input$password)
      if (is.null(u)) showModal(login_dialog("Wrong user or password")) else signed_in(u)
    })

    observeEvent(input$create_admin, {
      if (!identical(input$password, input$password2)) return(showModal(login_dialog("The passwords don't match")))
      result <- tryCatch(create_user(input$username, input$password, input$full_name, is_admin = TRUE),
                         error = function(e) e)
      if (inherits(result, "error")) return(showModal(login_dialog(conditionMessage(result))))
      signed_in(authenticate(input$username, input$password))
    })

    output$user_menu <- renderUI({
      u <- user()
      req(u)
      tagList(
        span(style = "color: white; margin-right: 10px;", icon("user"), u$username),
        actionLink(ns("change_password"), "Password", style = "color: white; margin-right: 10px;"),
        actionLink(ns("logout"), "Sign out", icon = icon("sign-out-alt"), style = "color: white;"))
    })

    # Signing out also makes this computer forget the user
    observeEvent(input$logout, {
      forget_user()
      session$reload()
    })

    observeEvent(input$change_password, {
      showModal(modalDialog(
        title = "Change password", size = "s", easyClose = TRUE,
        passwordInput(ns("old_password"), "Current password"),
        passwordInput(ns("new_password"), "New password (8 characters or more)"),
        passwordInput(ns("new_password2"), "Repeat the new password"),
        footer = tagList(modalButton("Cancel"), actionButton(ns("save_password"), "Save", class = "btn-success"))))
    })

    observeEvent(input$save_password, {
      u <- user()
      if (is.null(authenticate(u$username, input$old_password))) return(showNotification("The current password is wrong", type = "error"))
      if (!identical(input$new_password, input$new_password2)) return(showNotification("The passwords don't match", type = "error"))
      result <- tryCatch(set_password(u$id, input$new_password), error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      # a new password signs out the other computers; this one keeps remembering the user
      if (remember_enabled() && file.exists(remember_file())) remember_user(u$id)
      removeModal()
      showNotification("Password changed", type = "message")
    })

    user
  })
}
