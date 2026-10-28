test_that("users sign in with their password; inactive users can't", {
  fresh_store()
  expect_equal(user_count(), 0)
  id <- create_user("ana", "secret123", "Ana")
  expect_equal(authenticate("ANA", "secret123")$id, id)
  expect_null(authenticate("ana", "wrong-password"))
  expect_error(create_user("ana", "secret123"))          # unique
  expect_error(create_user("luis", "short"))             # 8+ characters
  update_user(id, active = FALSE)
  expect_null(authenticate("ana", "secret123"))
  update_user(id, active = TRUE)
  set_password(id, "another-password")
  expect_false(is.null(authenticate("ana", "another-password")))
  expect_false(grepl("another-password", paste(readLines(file.path(shared_dir, "accounts.json")), collapse = "")))
})

test_that("comments are visible to the author and to the groups they are shared with", {
  fresh_store()
  ana <- create_user("ana", "secret123"); luis <- create_user("luis", "secret123")
  eva <- create_user("eva", "secret123"); pedro <- create_user("pedro", "secret123")
  retina <- create_group("Retina"); cardio <- create_group("Cardio")
  set_user_groups(ana, c(retina, cardio)); set_user_groups(luis, retina); set_user_groups(eva, cardio)
  k <- variant_key("chr1", 100, "A", "G")
  add_comment(ana, k, "G", "Pathogenic", "to all my groups")
  add_comment(ana, k, "G", NA, "only retina", visibility = "group", group_id = retina)
  add_comment(ana, k, "G", NA, "private", visibility = "private")

  seen <- function(u) sort(visible_comments(u, k)$comment)
  expect_equal(seen(ana), sort(c("to all my groups", "only retina", "private")))
  expect_equal(seen(luis), sort(c("to all my groups", "only retina")))
  expect_equal(seen(eva), "to all my groups")
  expect_equal(seen(pedro), character(0))

  # group membership is checked when reading: leaving the group hides the comments
  set_user_groups(luis, integer(0))
  expect_equal(seen(luis), character(0))
})

test_that("each user writes only their own files in the shared folder", {
  fresh_store()
  ana <- create_user("ana", "secret123"); luis <- create_user("luis", "secret123")
  add_comment(ana, "k1", comment = "a"); add_comment(luis, "k1", comment = "b"); add_bookmark(ana, "s", "x")
  expect_setequal(list.files(file.path(shared_dir, "comments")), c("ana.jsonl", "luis.jsonl"))
  expect_equal(list.files(file.path(shared_dir, "sessions")), "ana.jsonl")
  # an incomplete last line (a write in progress in another computer) is ignored
  cat('{"id": "x", "variant_key": "k1"', file = file.path(shared_dir, "comments", "luis.jsonl"), append = TRUE)
  expect_equal(nrow(visible_comments(luis, "k1")), 1)
})

test_that("only the author deletes a comment; the newest classification wins", {
  fresh_store()
  ana <- create_user("ana", "secret123"); luis <- create_user("luis", "secret123")
  g <- create_group("Retina"); set_user_groups(ana, g); set_user_groups(luis, g)
  k <- variant_key("chr1", 100, "A", "G")
  add_comment(ana, k, classification = "Uncertain significance")
  Sys.sleep(1.1)
  add_comment(luis, k, classification = "Likely pathogenic", comment = "segregates")
  add_comment(luis, k, comment = "just a note")
  ann <- comment_annotations(ana)
  expect_equal(ann$classification, "Likely pathogenic")
  expect_equal(ann$n_comments, 3)
  id <- visible_comments(ana, k)$id[visible_comments(ana, k)$comment == "segregates"]
  expect_false(delete_comment(id, ana))
  expect_true(delete_comment(id, luis))
  expect_equal(comment_annotations(ana)$classification, "Uncertain significance")
  expect_error(add_comment(ana, k))   # empty
})

test_that("saved sessions belong to their user", {
  fresh_store()
  ana <- create_user("ana", "secret123"); luis <- create_user("luis", "secret123")
  add_bookmark(ana, "family 12", "1a2b")
  expect_equal(list_bookmarks(ana)$state_id, "1a2b")
  expect_equal(nrow(list_bookmarks(luis)), 0)
})

test_that("classifications of PriorR 3 can be imported", {
  fresh_store()
  ana <- create_user("ana", "secret123")
  f <- tempfile()
  writeLines(c("CHROM\tPOS\tREF\tALT\tClassification\tComments\tUser",
               "chr1\t100\tA\tG\tPathogenic\tseen before\tpablo"), f)
  expect_equal(import_legacy_classifications(f, ana), 1)
  cm <- visible_comments(ana, variant_key("chr1", 100, "A", "G"))
  expect_equal(cm$classification, "Pathogenic")
  expect_match(cm$comment, "PriorR 3 classification by pablo")
})
