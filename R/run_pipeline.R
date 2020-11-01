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

  # initialize
  experiment_name <- basename(data_path)

  create_data_subdir(data_path)

  config <- load_config(data_path)

  # start analysis
  cs <- load_fcs(data_path)

  bead_index <- grepl(config$beads_pattern, flowCore::sampleNames(cs))

  if (sum(bead_index) == 1) {

    cs_beads <- cs[bead_index]
    cs <- cs[!bead_index]

    t_fun <- generate_transformation(cs_beads)

  } else if (sum(bead_index) > 1) {

    stop("More than one bead sample found in the dataset. Please restricit the
         bead_pattern to one sample.")
  }

  # gating
  gt_file <- system.file("tools", "gt_samples.csv", package = "expressalyzr",
                         mustWork = TRUE)

  gt <- openCyto::gatingTemplate(gt_file)
  gs <- flowWorkspace::GatingSet(cs)

  openCyto::register_plugins(fun = mixture_gate, methodName = "mixture_gate",
                             dep = c("Rmixmod"), "gating")

  openCyto::gt_gating(gt, gs)

  data_dt <- flowWorkspace::gs_pop_get_data(gs, y = "singlets")

  trans_dt <- apply_transform(data_dt)

}
