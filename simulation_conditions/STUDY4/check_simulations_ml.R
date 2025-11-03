# check_simulations_ml.R — quick per-unit diagnostics
`%||%` <- function(a, b) if (!is.null(a)) a else b

compute_residuals_TRUE <- function(y, mu, Phi) {
  Tn <- nrow(y)
  if (Tn < 2) {
    return(list(res1 = numeric(0), res2 = numeric(0)))
  }
  muM <- matrix(mu, nrow = Tn - 1, ncol = 2, byrow = TRUE)
  pred <- muM + y[1:(Tn - 1), ] %*% t(Phi)
  res <- y[2:Tn, ] - pred
  list(res1 = res[, 1], res2 = res[, 2])
}

plot_unit <- function(df, mu_i, Phi, skew_dir, cond_lbl, rep_lbl, unit_lbl) {
  y <- as.matrix(df[, c("y1", "y2")])
  res <- compute_residuals_TRUE(y, mu_i, Phi)

  old <- par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old))
  plot(df$t, df$y1, type = "l", main = "Series y1", xlab = "t", ylab = "y1")
  plot(df$t, df$y2, type = "l", main = "Series y2", xlab = "t", ylab = "y2")
  plot(df$y1, df$y2, pch = 19, col = "#00000044", main = "Scatter y1 vs y2", xlab = "y1", ylab = "y2")

  if (length(res$res1) > 0) {
    hist(res$res1, breaks = "fd", main = paste0("Hist res1 (", ifelse(skew_dir[1] > 0, "Right", "Left"), ")"), xlab = "res1")
    hist(res$res2, breaks = "fd", main = paste0("Hist res2 (", ifelse(skew_dir[2] > 0, "Right", "Left"), ")"), xlab = "res2")
    qqnorm(res$res1, main = "QQ-Norm y1")
    qqline(res$res1)
    qqnorm(res$res2, main = "QQ-Norm y2")
    qqline(res$res2)
    acf(res$res1, main = "ACF res1", na.action = na.pass)
    acf(res$res2, main = "ACF res2", na.action = na.pass)
  }
  mtext(sprintf("%s | Rep %d | Unit %d", cond_lbl, rep_lbl, unit_lbl), outer = TRUE, line = 0, cex = 1.1, font = 2)
}

run_post_sim_checks_var1_ml <- function(data_dir, checks_dir, n_units_to_plot = 3) {
  files <- list.files(data_dir, "^sim_dataML_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(files)) {
    message("No ML simulation files found.")
    return(invisible())
  }

  cids <- sort(unique(sub("^.*_cond(\\d+)_rep.*", "\\1", basename(files))))
  dir.create(checks_dir, FALSE, TRUE)

  for (cid in cids) {
    reps <- list.files(data_dir, sprintf("^sim_dataML_cond%03d_rep\\d+\\.rds$", as.integer(cid)), full.names = TRUE)
    if (!length(reps)) next
    set.seed(as.integer(cid))
    reps_to_check <- if (length(reps) > 3) sample(reps, 3) else reps

    pdf(file.path(checks_dir, sprintf("checks_ml_cond_%03d.pdf", as.integer(cid))), width = 11, height = 8.5)
    for (rf in reps_to_check) {
      dat <- try(readRDS(rf), silent = TRUE)
      if (inherits(dat, "try-error")) next
      Phi <- dat$phi_matrix
      mu_mat <- dat$true_params$mu_mat
      skew_dir <- dat$true_params$skew_dir %||% c(1, 1)
      cond_lbl <- sprintf("Cond %03d (ML-Exponential)", dat$condition_id)
      # choose a few units
      units <- sort(sample(unique(dat$data$i), min(n_units_to_plot, dat$N)))
      for (u in units) {
        dfu <- dat$data[dat$data$i == u, c("t", "y1", "y2")]
        plot_unit(dfu, mu_mat[u, ], Phi, skew_dir, cond_lbl, dat$rep_i, u)
      }
    }
    dev.off()
  }
  message("ML simulation checks complete.")
}
