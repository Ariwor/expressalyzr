#' Run the data processing pipline
#'
#' Run the data processing pipline that links together other functions in the
#'  expressalyzr package.
#'
#' @return
#' @export
#'
run_pipeline <- function(data_path, view_config = TRUE) {

  # initialize
  experiment_name <- basename(data_path)

  create_data_subdir(data_path)

  config_file_path <- file.path(data_path, "config.yml")
  config <- load_config(config_file_path, view_config)

  # start analysis
  cs <- load_fcs(data_path)

  bead_index <- grepl(config$beads_pattern, flowCore::sampleNames(cs))

  if (sum(bead_index) == 1) {

    cs_beads <- cs[bead_index]
    cs <- cs[!bead_index]

    if (config$mefl_transform) {
      t_fun <- generate_transformation(cs_beads)
      trans <- TRUE
    } else {
      trans <- FALSE
    }
  } else if (sum(bead_index) > 1) {

    stop("More than one bead sample found in the dataset. Please restricit the
         bead_pattern to one sample.")

  } else {

    message("Bead sample not found. Values will not be transformed to MEFL.")
    trans <- FALSE
  }

  # gating
  gt_file <- system.file("tools", "gt_samples.csv", package = "expressalyzr",
                         mustWork = TRUE)

  gt <- openCyto::gatingTemplate(gt_file)
  gs <- flowWorkspace::GatingSet(cs)

  openCyto::register_plugins(fun = mixture_gate, methodName = "mixture_gate",
                             dep = c("Rmixmod"), "gating")

  openCyto::gt_gating(gt, gs)

  data_cs <- flowWorkspace::gs_pop_get_data(gs, y = "singlets")

  if (trans) {
    data_cs <- apply_transform(data_cs, t_fun)
  }

  # compensation
  cont_index <- grepl(config$controls_pattern, flowCore::sampleNames(data_cs))

  if (sum(cont_index) > 2) {

    cont_cs <- data_cs[cont_index]

    so_mat_path <- file.path(data_path, "so_mat.RData")

    if (!file.exists(so_mat_path) || config$redo_comp) {

      so_mat <- spillover_matrix(cont_cs,
                                 config$controls_index,
                                 config$channel_pattern,
                                 config$density_th,
                                 config$manual_comp)

      save(so_mat, file = so_mat_path)
    } else (
      load(so_mat_path)
    )

    data_cs <- flowWorkspace::compensate(data_cs, so_mat)

    cont_cs <- data_cs[cont_index]
    data_cs <- data_cs[!cont_index]

    cont <- TRUE
  } else {
    cont <- FALSE
    message("No control samples found. Proceeding with out compensation.")
  }
  # convert cytoset to data table
  data_dt <- cs_to_dt(data_cs)

  chs <- colnames(data_dt)[grepl("FL", colnames(data_dt))]

  # gate populations
  data_dt[, no_negative := rowSums(data_dt[, chs, with = FALSE] <= 0) == 0]

  if (cont) {
    cont_dt <- cs_to_dt(cont_cs)
    neg_dt <- cs_to_dt(cont_cs[config$controls_index[1]])

    if (config$manual_cutoff) {
      select_chs <- chs[grepl(config$channel_pattern, chs)]
      cont_sub_dt <- cont_dt[, c("File", select_chs), with = FALSE]
      config$bg_cutoff <- adjust_threshold(cont_sub_dt, config$bg_cutoff)
      write_value <- paste("bg_cutoff:", config$bg_cutoff)
      write_config(config_file_path, "default", write_value)
    }

    neg_dt[, lapply(mget(chs), quantile, config$bg_cutoff),
           by = .(File)]

    bg_dt <- neg_dt[, lapply(mget(chs), quantile, config$bg_cutoff)]

    f_pos <- function(channel, values) values >= bg_dt[[channel]]

    pos_chs <- paste0(chs, "_pos")
    data_dt[, (pos_chs) := lapply(chs, function(ch) f_pos(ch, get(ch))), by = .(File)]
  }

  # assign experimental specifications
  s_file_path <- file.path(data_path, config$spec_file)

  if (file.exists(s_file_path)) {
    s_file <- data.table::fread(s_file_path)
    data_dt[, ID := gsub("^0(\\d{1})-.*-([A-Z]{1}\\d{1,2})\\.fcs$", "\\1-\\2", File)]
    data_dt <- merge(data_dt, s_file, by = config$merge_by, all.x = TRUE)
  }

  data.table::fwrite(data_dt, file = file.path(data_path, paste0(experiment_name, ".csv")))
  return(data_dt)
}
