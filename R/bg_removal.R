#'
#'
assign_bg <- function(x, n_comp = NULL, log = TRUE, inspect = FALSE, rm = 1) {
  if (log) {
    x <- log(x)
  }
  cv_x <- sd(x) / mean(x) * 100

  if (is.null(n_comp)) {
    if (cv_x > 18) {
      n_comp <- 3
    } else {
      n_comp <- 2
    }
  }

  if (n_comp == 2) {
    fast <- TRUE
  } else {
    fast <- FALSE
  }
  fit <- mixtools::normalmixEM(x, k = n_comp, maxit = 50000, epsilon = 1e-08, fast = fast)
  bg_index <- which(fit$mu %in% head(sort(fit$mu), rm))
  if (inspect) {
    mixtools::plot.mixEM(fit, whichplots = 2, n = 120)
    title(sub = paste("\t\t\t\t\t\t\t\t\t\t\t\tCV = ", cv_x,
                      "\n\t\t\t\t\t\t\t\t\t\t\t\tAIC = ", 2 * fit$loglik + (n_comp * 3 - 1) * 2,
                      "\n\t\t\t\t\t\t\t\t\t\t\t\tBIC = ", 2 * fit$loglik + (n_comp * 3 - 1) * log(length(x))),
          xlab = "")
  }

  if (length(bg_index) < 2) {
    bg <- fit$posterior[, bg_index] > runif(length(x))
  } else {
    bg <- apply(fit$posterior[, bg_index], 1, sum) > runif(length(x))
  }
  return(bg)
}
