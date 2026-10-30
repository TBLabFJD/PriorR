
test_that("first run creates the administrator, then users sign in", {
  fresh_store()
  testServer(auth_server, {
    session$setInputs(username = "admin", full_name = "Admin", password = "secret123", password2 = "secret123")
    session$setInputs(create_admin = 1)
    expect_equal(user()$username, "admin")
    expect_true(user()$is_admin)
  })
  testServer(auth_server, {
    session$setInputs(username = "admin", password = "wrong-password", login = 1)
    expect_null(user())
    session$setInputs(username = "admin", password = "secret123", login = 2)
    expect_equal(user()$username, "admin")
  })
})

test_that("loading: nothing is loaded without a user", {
  fresh_store()
  testServer(load_server, args = list(user = reactiveVal(NULL)), {
    session$setInputs(file = open_file_value(make_variants(3)), profile = "builtin:None", run = 1)
    expect_false(state$running)
  })
})

wait_loaded <- function(session, loaded) {
  for (i in 1:60) {
    if (!is.null(loaded())) return(loaded())
    Sys.sleep(0.5)
    session$elapse(700)
  }
  stop("not loaded")
}

test_that("loading: the same way for any file, with the load filters of the profile", {
  fresh_store()
  create_user("ana", "secret123")
  me <- reactiveVal(authenticate("ana", "secret123"))
  df <- make_variants(6, Genomic_region = c("EXONIC", "INTRONIC", "EXONIC", "EXONIC", "UPSTREAM", "SPLICING"),
                      gnomADg_AF_popmax = c(0.001, 0.001, 0.5, NA, 0.001, 0.001))
  testServer(load_server, args = list(user = me), {
    session$setInputs(file = open_file_value(df), profile = "builtin:None", regions = NULL, canonical = FALSE,
                      af_gnomadg_popmax = NA, af_gnomadg = NA, af_gnomade = NA, af_fjd = NA, min_cadd = NA,
                      excl_benign = FALSE, panels = NULL, rescue_clinvar = FALSE, rescue_spliceai = NA, max_rows = 1000)
    session$setInputs(run = 1)
    l <- wait_loaded(session, loaded)
    expect_equal(nrow(l$df), 6)                      # None: everything
    expect_equal(l$profile$name, "None")
    loaded(NULL)
    session$setInputs(profile = "builtin:WGS", regions = c("EXONIC", "SPLICING"), af_gnomadg_popmax = 0.01)
    session$setInputs(run = 2)
    l <- wait_loaded(session, loaded)
    expect_equal(l$df$POS, df$POS[c(1, 4, 6)])       # coding/splicing and rare
    expect_equal(l$profile$name, "WGS")
    # genes typed and of a saved panel, without the excluded panel
    save_custom_panel("Excl", "GENE9")
    loaded(NULL)
    session$setInputs(regions = NULL, af_gnomadg_popmax = NA, genes = "GENE1, GENE9", excluded_panels = "custom:Excl")
    session$setInputs(run = 3)
    l <- wait_loaded(session, loaded)
    expect_equal(nrow(l$df), 6)                      # all are GENE1
    expect_equal(unlist(settings()$genes), c("GENE1", "GENE9"))
  })
})

test_that("SNV module: a loaded file gets the filters of its profile; filter, classify and annotate", {
  fresh_store()
  create_user("ana", "secret123")
  me <- reactiveVal(authenticate("ana", "secret123"))
  df <- make_variants(6, CADD_PHRED = c(30, 5, 20, NA, 16, 40), CLNSIG = c(NA, NA, "Benign", NA, "Pathogenic", NA),
                      Genomic_region = c("EXONIC", "EXONIC", "EXONIC", "INTRONIC", "SPLICING", NA))
  testServer(snv_server, args = list(user = me), {
    session$setInputs(consequence = c("missense_variant", "NA"), region = c(genomic_regions, "NA"),
                      type = "SNV", canonical = FALSE, exclude_predicted = FALSE, genes = "",
                      max_gnomadg_af = 1, max_gnomadg_popmax = 1, max_gnomade_af = 1, max_fjd_af = 1,
                      min_cadd = "0", exclude_benign = FALSE, inheritance = "", glow = -1, orphanet = NULL,
                      panels = NULL, amber = FALSE, excluded_panels = NULL, phenotype = NULL, proband = "",
                      other_test_filter = FALSE, view = "filters", color_priority = FALSE, annotate = FALSE, columns = NULL)
    loader$loaded(list(df = df, profile = builtin_profiles()$None, restoring = FALSE))
    session$flushReact()
    expect_equal(nrow(filtered()), 6)                 # default: no filters
    expect_equal(applied_profile(), "None")

    # the WES profile as the sidebar sets it (CADD >= 15 as a number, benign out, coding/splicing)
    session$setInputs(region = c("EXONIC", "SPLICING"), min_cadd = "15", exclude_benign = TRUE, view = "priority",
                      columns = c("SYMBOL", "CADD_PHRED"))
    out <- filtered()
    expect_equal(out$CADD_PHRED, c(16, 30))
    expect_equal(attr(out, "display_columns"), c("SYMBOL", "CADD_PHRED"))

    # the current filters and columns are saved as a profile of the user
    session$setInputs(profile_name = "Mi retina", save_profile = 1)
    saved <- user_profiles("ana")[["Mi retina"]]
    expect_equal(unlist(saved$filters$region), c("EXONIC", "SPLICING"))
    expect_equal(saved$filters$min_cadd, 15)
    expect_equal(unlist(saved$columns), c("SYMBOL", "CADD_PHRED"))
    expect_null(saved$filters$type)                   # all types selected = no filter

    # Classify row 2 -> the comment is saved for that variant
    session$setInputs(classify = 2)
    info$open(list(key = variant_key(out$CHROM[2], out$POS[2], out$REF[2], out$ALT[2]), gene = "GENE1", label = "x"))
    session$setInputs(`info-comments-classification` = "Likely pathogenic", `info-comments-comment` = "check",
                      `info-comments-share` = "groups", `info-comments-save` = 1)
    session$setInputs(annotate = TRUE)
    out <- filtered()
    expect_equal(out$Classification, c(NA, "Likely pathogenic"))
    expect_equal(out$Comments, c(0L, 1L))
  })
})

test_that("gene panels of the library filter the table", {
  fresh_store()
  save_custom_panel("Mi panel", c("GENE2", "GENE3"))
  me <- reactiveVal({ create_user("ana", "secret123"); authenticate("ana", "secret123") })
  df <- make_variants(4, SYMBOL = c("GENE1", "GENE2", "GENE3", "GENE4"))
  testServer(snv_server, args = list(user = me), {
    session$setInputs(consequence = NULL, region = NULL, type = NULL, canonical = FALSE, exclude_predicted = FALSE,
                      genes = "", max_gnomadg_af = 1, max_gnomadg_popmax = 1, max_gnomade_af = 1, max_fjd_af = 1,
                      min_cadd = "0", exclude_benign = FALSE, inheritance = "", glow = -1, orphanet = NULL,
                      panels = "custom:Mi panel", amber = FALSE, excluded_panels = NULL, phenotype = NULL, proband = "",
                      other_test_filter = FALSE, view = "filters", color_priority = FALSE, annotate = FALSE, columns = NULL)
    loader$loaded(list(df = df, profile = builtin_profiles()$None, restoring = TRUE))
    session$flushReact()
    expect_equal(filtered()$SYMBOL, c("GENE2", "GENE3"))
  })
})

test_that("only administrators manage users", {
  fresh_store()
  admin <- create_user("admin", "secret123", is_admin = TRUE)
  testServer(admin_server, args = list(user = reactiveVal(authenticate("admin", "secret123"))), {
    session$setInputs(new_group = "Retina", new_group_description = "", create_group = 1)
    session$setInputs(new_username = "luis", new_full_name = "Luis", new_password = "secret123",
                      new_is_admin = FALSE, new_groups = list_groups()$id, create_user = 1)
    expect_equal(user_groups(authenticate("luis", "secret123")$id)$name, "Retina")
  })
  luis <- authenticate("luis", "secret123")
  testServer(admin_server, args = list(user = reactiveVal(luis)), {
    session$setInputs(new_group = "Hacked", new_group_description = "", create_group = 1)
    expect_false("Hacked" %in% list_groups()$name)
  })
})

test_that("remember me: a computer that remembers its user goes straight in until signing out", {
  fresh_store()
  withr::local_envvar(PRIORR_LOCAL = "1", PRIORR_LOGIN_DIR = tempfile("login-"))
  create_user("ana", "secret123")
  testServer(auth_server, {
    expect_null(user())                                   # first time: the sign in dialog
    session$setInputs(username = "ana", password = "secret123", remember = TRUE, login = 1)
    expect_equal(user()$username, "ana")
  })
  expect_true(file.exists(remember_file()))
  expect_false(grepl(jsonlite::fromJSON(remember_file())$token, paste(readLines(file.path(shared_dir, "accounts.json")), collapse = "")))
  testServer(auth_server, expect_equal(user()$username, "ana"))     # next time: straight in
  set_password("ana", "new-password")                               # e.g. reset by an administrator
  testServer(auth_server, expect_null(user()))                      # every computer forgets the user
})

test_that("remember me: signing out forgets the user; disabled when not running locally", {
  fresh_store()
  withr::local_envvar(PRIORR_LOCAL = "1", PRIORR_LOGIN_DIR = tempfile("login-"))
  create_user("ana", "secret123")
  remember_user("ana")
  expect_equal(remembered_user()$username, "ana")
  forget_user()
  expect_null(remembered_user())
  remember_user("ana")
  withr::local_envvar(PRIORR_LOCAL = "")
  expect_null(remembered_user())
})

test_that("profiles in the load dialog: save, upload a colleague's profile, delete", {
  fresh_store()
  me <- reactiveVal({ create_user("ana", "secret123"); authenticate("ana", "secret123") })
  view <- function() list(filters = list(min_cadd = 20, region = list("EXONIC")), columns = list("SYMBOL", "HGVSc"))
  testServer(load_server, args = list(user = me, current_view = view), {
    session$setInputs(profile = "builtin:WES", regions = c("EXONIC", "SPLICING"), canonical = TRUE,
                      af_gnomadg_popmax = 0.005, af_gnomadg = NA, af_gnomade = NA, af_fjd = NA, min_cadd = 10,
                      excl_benign = TRUE, rescue_clinvar = TRUE, rescue_spliceai = 0.5, max_rows = 5000,
                      genes = "", panels = "custom:none", amber = FALSE, excluded_panels = NULL)
    session$setInputs(profile_name = "Mi exoma", save_profile = 1)
    saved <- user_profiles("ana")[["Mi exoma"]]
    expect_equal(saved$load$max_af$gnomADg_AF_popmax, 0.005)      # the load filters of the dialog
    expect_equal(saved$load$min_cadd, 10)
    expect_equal(saved$filters$min_cadd, 20)                      # the sidebar filters
    expect_equal(unlist(saved$columns), c("SYMBOL", "HGVSc"))
    expect_true("user:Mi exoma" %in% profile_choices("ana"))

    # a profile downloaded by a colleague
    f <- tempfile(fileext = ".json")
    jsonlite::write_json(list(name = "Cardio de Luis", load = list(regions = list("EXONIC")), filters = list(min_cadd = 25)),
                         f, auto_unbox = TRUE)
    session$setInputs(profile_file = data.frame(name = "cardio.json", datapath = f))
    expect_equal(user_profiles("ana")[["Cardio de Luis"]]$filters$min_cadd, 25)
    expect_equal(user_profiles("ana")[["Cardio de Luis"]]$filters$max_gnomadg_af, 1)

    session$setInputs(profile = "user:Mi exoma", delete_profile = 1)
    expect_false("Mi exoma" %in% names(user_profiles("ana")))
  })
})
