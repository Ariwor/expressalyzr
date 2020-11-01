#### errors if non-existant file path or no .fcs files in the directory ####

test_that("non-existing file paths throw an error", {
  test_dir_path <- file.path(".", "does_not_exist")
  # create dummy directory just to remove again
  expect_true(dir.create(test_dir_path))
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
  expect_error(load_fcs(test_dir_path), "^.* not exist .*$")
})

test_that("existing file paths with no .fcs files throw an error", {
  test_dir_path <- file.path(".", "no_fcs_files")
  # create dummy directory for testing
  expect_true(dir.create(test_dir_path))
  expect_true(dir.exists(test_dir_path))
  expect_error(load_fcs(test_dir_path), "^.* not contain any .*$")
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
})

#### return value  ####

test_that("a cytoset is returned", {
  test_file_path <- system.file("extdata", "example_fcs_files",
                                package = "expressalyzr", mustWork = TRUE)
  expect_equal(is(load_fcs(test_file_path)), c("cytoset", "flowSet"))
})
