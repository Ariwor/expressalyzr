#'
mixture_gate <- function(fr, pp_res, channels = NA, filterId = "",
                         n_samples = NULL, alg = "CEM", th = 2e-14,
                         n_clusters = 1:5, log_t = FALSE) {
  print(th)
  max_r <- range(fr)[2, ]
  dat <- as.data.frame(flowCore::exprs(fr))
  dat <- dat[dat$`FSC-A` <= max_r[[1]] & dat$`SSC-A` <= max_r[[2]], ]

  if (!is.null(n_samples) && nrow(dat) > n_samples) {
    sample_ind <- sample.int(nrow(dat), n_samples)
    dat <- dat[sample_ind, ]
  }

  if (log_t) {
    dat <- log(dat[apply(dat > 0, 1, all), ])
  }

  mods <- Rmixmod::mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Ck")

  strat <- Rmixmod::mixmodStrategy(algo = alg,
                                   nbTry = 50,
                                   initMethod = "smallEM",
                                   epsilonInAlgo = 1e-06)

  fit <- Rmixmod::mixmodCluster(dat,
                                dataType = "quantitative",
                                nbCluster = n_clusters,
                                models = mods,
                                strategy = strat,
                                criterion = "BIC")

  clu_i <- which.max(fit["parameters"]["proportions"])

  mu <- fit["parameters"]["mean"][clu_i, ]
  Sigma <- fit["parameters"]["variance"][[clu_i]]

  clu_s <- dat[fit["partition"] == clu_i, ]

  xs <- seq(min(dat[, 1]), max(dat[, 1]), length.out = 250)
  ys <- seq(min(dat[, 2]), max(dat[, 2]), length.out = 250)

  grid <- as.matrix(expand.grid(xs, ys))

  prob <- Rfast::dmvnorm(grid, mu, Sigma)
  gated <- grid[prob > th,]
  gate <- gated[chull(gated), ]

  if (log_t) {
    gate <- exp(gate)
  }

  colnames(gate) <- channels

  return(flowCore::polygonGate(.gate = gate, filterId = filterId))
}
