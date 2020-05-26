context("error handling")

test_that("non-existing file paths throw an error", {
  test_file_path <- "./doesnotexist/"
  expect_error(load_fcs(test_file_path), "^.* not exist .*$")
})
