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
  fcs_files <- list.files(path = file_path, pattern = "^.*\\.fcs$")
  if (identical(fcs_files, character(0))) {
    stop("The specified file path does either not exist or not contain any .fcs files.")
  } else {
    gating_set <- ncdfFlow::read.ncdfFlowSet(fcs_files)
  }
  return(gating_set)
}
