# Run all the tests from the repository root:  Rscript tests/run_tests.R
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
