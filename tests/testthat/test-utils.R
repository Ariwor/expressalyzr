#### test create_data_subdir function ####

test_that("a new directory with the name 'data' is created", {
  # setup
  test_dir_path <- file.path(".", "dir_without_data_subdir")
  dummy_file <- file.path(test_dir_path, "dummy.fcs")
  dummy_subdir <- file.path(test_dir_path, "dummy_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(file.create(dummy_file))
  expect_true(dir.create(dummy_subdir))
  # test
  expect_true(create_data_subdir(test_dir_path))
  expect_true(dir.exists(file.path(test_dir_path, "data")))
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})

test_that("files and folders are moved from the old to the new data directory", {
  # setup
  test_dir_path <- file.path(".", "dir_for_moving_data")
  dummy_file <- file.path(test_dir_path, "dummy.fcs")
  dummy_subdir <- file.path(test_dir_path, "dummy_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(file.create(dummy_file))
  expect_true(dir.create(dummy_subdir))
  # test
  expect_true(create_data_subdir(test_dir_path))
  expect_true(file.exists(file.path(test_dir_path, "data", "dummy.fcs")))
  expect_true(dir.exists(file.path(test_dir_path, "data", "dummy_subdir")))
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})

test_that("nothing happens if data subdirectory already exists", {
  # setup
  test_dir_path <- file.path(".", "dir_with_data_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(dir.create(file.path(test_dir_path, "data")))
  # test
  expect_null(create_data_subdir(test_dir_path))
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})
