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
