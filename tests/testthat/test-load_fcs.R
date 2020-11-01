context("test fcs file loading")

#### errors if non-existant file path or no .fcs files in the directory ####

test_that("non-existing file paths throw an error", {
  test_file_path <- "./does_not_exist/"
  # create dummy directory just to remove again
  expect_equal(dir.create(test_file_path), TRUE)
  expect_equal(unlink(test_file_path, recursive = TRUE), 0)
  expect_equal(dir.exists(test_file_path), FALSE)
  expect_error(load_fcs(test_file_path), "^.* not exist .*$")
})

test_that("existing file paths with no .fcs files throw an error", {
  test_file_path <- "./no_fcs_files"
  # create dummy directory for testing
  expect_equal(dir.create(test_file_path), TRUE)
  expect_equal(dir.exists(test_file_path), TRUE)
  expect_error(load_fcs(test_file_path), "^.* not contain any .*$")
  # clean up
  expect_equal(unlink(test_file_path, recursive = TRUE), 0)
})

#### return value  ####

test_that("a gating set is returned", {
  test_file_path <- system.file("extdata", "example_fcs_files", package = "expressalyzr", mustWork = TRUE)

})
