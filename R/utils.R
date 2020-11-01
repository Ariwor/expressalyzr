#' Load .fcs files
#'
#' Load .fcs files into GatingSet object using read.ncdfFlowSet function.
#'
#' @param file_path The path to the .fcs files that are to be loaded.
#'
#' @return A GatingSet object.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_fcs_files", package = "expressalyzr", mustWork = TRUE)
#' gating_set <- load_fcs(path)
load_fcs <- function(file_path) {
  fcs_files <- list.files(path = file_path, pattern = "^.*\\.fcs$", full.names = TRUE)
  if (identical(fcs_files, character(0))) {
    stop("The specified file path does either not exist or does not contain any .fcs files.")
  } else {
    gating_set <- flowWorkspace::load_cytoset_from_fcs(files = fcs_files)
  }
  return(gating_set)
}

# this function causes a segfault form C stack overflow on my computer
load_fcs_ncdf <- function(file_path) {
  fcs_files <- list.files(path = file_path, pattern = "^.*\\.fcs$", full.names = TRUE)
  if (identical(fcs_files, character(0))) {
    stop("The specified file path does either not exist or does not contain any .fcs files.")
  } else {
    gating_set <- ncdfFlow::read.ncdfFlowSet(fcs_files)
  }
  return(gating_set)
}

#' Create data subdirectory
#'
#' This is a utility function for creating a subdirectory for the data in the
#' data path and moving the data files into the new directory. This way the
#' analysis results and data stay separated.
#'
#' @param data_path The path to the original data directory.
#'
#' @return
#' @export
#'
#' @example
#'
create_data_subdir <- function(data_path) {

  new_data_path <- file.path(data_path, "data")

  if (!dir.exists(new_data_path)) {
    # create new folder
    data_files <- list.files(data_path)
    create_out <- dir.create(new_data_path)

    # move files
    data_files_from <- file.path(data_path, data_files)
    data_files_to <- file.path(new_data_path, data_files)

    rename_out <- file.rename(data_files_from, data_files_to)

    output <- create_out & all(rename_out)
  } else {
    output <- NULL
  }

  return(output)
}

#' Create and/or load configuration file
#'
#' A configuration file is created from a template and opened for editing.
#' After editing is complete or if a configuration file already exists
#' the file is loaded.
#'
#' @param
#' @export
#'
#' @example
#'
manage_config <- function(data_path) {

  config_file_path <- file.path(data_path, "config.yml")

  if (!exists(config_file_path)) {

  }
}
