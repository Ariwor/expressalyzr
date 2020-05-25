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
#' gating_set <- load_fcs(<path to fcs files>)
load_fcs <- function(file_path) {
  fcs_files <- list.files(path = file_path, pattern = "^.*\\.fcs$")
  gating_set <- ncdfFlow::read.ncdfFlowSet(fcs_files)
}
