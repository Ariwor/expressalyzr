#' Run the data processing pipline
#'
#' Run the data processing pipline that links together other functions in the
#'  expressalyzr package.
#'
#' @return
#' @export
#'
#' @examples
#' run_pipeline()
run_pipeline <- function(data_path) {

  experiment_name <- basename(data_path)

  create_data_subdir(data_path)



  cs <- load_fcs(data_path)
}
