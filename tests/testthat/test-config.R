test_that("priorr.cfg: key: value lines, comments and Windows paths", {
  f <- tempfile()
  writeLines(c("# comment", "", "shared_folder: Z:\\priorr_datos", "data_folders: Genomes=/data/g; Home=~", "threads: 8"), f)
  cfg <- read_config(f)
  expect_equal(cfg$shared_folder, "Z:\\priorr_datos")
  expect_equal(cfg$threads, "8")
  roots <- parse_roots(cfg$data_folders)
  expect_equal(names(roots), c("Genomes", "Home"))
  expect_equal(unname(roots["Home"]), path.expand("~"))
  expect_equal(read_config(file.path(tempdir(), "missing.cfg")), list())
})

test_that("settings: environment variable > computer config > lab config > default", {
  user_cfg <- tempfile(); lab_home <- tempfile(); dir.create(lab_home)
  withr::local_envvar(PRIORR_USER_CONFIG = user_cfg, PRIORR_WGS_THREADS = "")
  old_home <- priorr_home; priorr_home <<- lab_home; on.exit(priorr_home <<- old_home)
  expect_equal(setting("threads", "PRIORR_WGS_THREADS", "4"), "4")
  writeLines("threads: 6", file.path(lab_home, "priorr.cfg"))
  expect_equal(setting("threads", "PRIORR_WGS_THREADS", "4"), "6")
  writeLines("threads: 2", user_cfg)
  expect_equal(setting("threads", "PRIORR_WGS_THREADS", "4"), "2")
  withr::local_envvar(PRIORR_WGS_THREADS = "12")
  expect_equal(setting("threads", "PRIORR_WGS_THREADS", "4"), "12")
})

test_that("first-run setup only on a user's computer without a configured shared folder", {
  lab_home <- tempfile(); dir.create(lab_home)
  old_home <- priorr_home; priorr_home <<- lab_home; on.exit(priorr_home <<- old_home)
  withr::local_envvar(PRIORR_USER_CONFIG = tempfile(), PRIORR_SHARED_DIR = "", PRIORR_LOCAL = "1")
  expect_false(shared_folder_configured())
  withr::local_envvar(PRIORR_LOCAL = "")
  expect_true(shared_folder_configured())                       # a server: no setup
  withr::local_envvar(PRIORR_LOCAL = "1")
  writeLines("shared_folder: /mnt/lab", file.path(lab_home, "priorr.cfg"))
  expect_true(shared_folder_configured())                       # configured for the whole lab
})

test_that("the setup saves the shared folder in the computer config", {
  lab_home <- tempfile(); dir.create(lab_home)
  old_home <- priorr_home; priorr_home <<- lab_home
  old_shared <- shared_dir
  on.exit({ priorr_home <<- old_home; shared_dir <<- old_shared })
  user_cfg <- tempfile()
  withr::local_envvar(PRIORR_USER_CONFIG = user_cfg, PRIORR_SHARED_DIR = "", PRIORR_LOCAL = "1",
                      PRIORR_LOCAL_SHARED = tempfile("local-shared-"))
  testServer(setup_server, {
    session$setInputs(mode = "local", save = 1)
  })
  expect_equal(read_config(user_cfg)$shared_folder, local_shared_dir())
  expect_equal(shared_dir, local_shared_dir())
  expect_true(shared_folder_configured())
})

test_that("folders of the lab configuration can be relative to the PriorR folder", {
  expect_true(is_absolute_path("/data/x")); expect_true(is_absolute_path("~/x"))
  expect_true(is_absolute_path("D:\\priorr")); expect_true(is_absolute_path("\\\\server\\genetics"))
  expect_false(is_absolute_path("shared")); expect_false(is_absolute_path("data/shared"))
  old <- priorr_home
  priorr_home <<- tempfile("priorr-home-"); dir.create(priorr_home)
  writeLines("shared_folder: shared", file.path(priorr_home, "priorr.cfg"))
  withr::with_envvar(c(PRIORR_SHARED_DIR = "", PRIORR_USER_CONFIG = tempfile()), {
    expect_equal(setting("shared_folder", "PRIORR_SHARED_DIR", "x"), file.path(priorr_home, "shared"))
  })
  priorr_home <<- old
})
