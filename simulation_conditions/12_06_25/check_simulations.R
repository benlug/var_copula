###########################################################################
# check_simulations.R
# Visual sanity checks – produces one PDF per condition (max 5 reps each)
###########################################################################

# Needed base packages only (graphics, grDevices, stats)

compute_residuals_TRUE <- function(y, mu, phi) {
  Tn <- nrow(y)
  if (Tn < 2) {
    return(list(res1 = numeric(0), res2 = numeric(0)))
  }
  res <- y[2:Tn, ] - (matrix(mu, nrow = Tn - 1, ncol = 2, byrow = TRUE) +
    y[1:(Tn - 1), ] %*% t(phi))
  list(res1 = res[, 1], res2 = res[, 2])
}

plot_diagnostics <- function(df, true_par, cond_lbl, rep_lbl) {
  ymat <- as.matrix(df[, c("y1", "y2")])
  res <- compute_residuals_TRUE(ymat, true_par$mu, true_par$phi)

  old <- par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
  on.exit(par(old))

  plot(df$t, df$y1, type = "l", main = "Series y1")
  plot(df$t, df$y2, type = "l", main = "Series y2")
  plot(df$y1, df$y2,
    pch = 19, col = "#00000044",
    main = "Scatter"
  )

  hist(res$res1, breaks = "fd", main = "Resid y1")
  hist(res$res2, breaks = "fd", main = "Resid y2")
  qqnorm(res$res1, main = "QQ y1")
  qqline(res$res1)
  qqnorm(res$res2, main = "QQ y2")
  qqline(res$res2)

  acf(res$res1, main = "ACF res1")
  acf(res$res2, main = "ACF res2")
  mtext(sprintf("%s  |  Rep %d", cond_lbl, rep_lbl),
    outer = TRUE, line = -1.2, cex = 1.1, font = 2
  )
}

run_post_sim_checks_var1 <- function(data_dir, checks_dir) {
  files <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(files)) {
    message("No simulation files found.")
    return(invisible())
  }

  cond_ids <- sort(unique(sub(
    "^.*_cond(\\d+)_rep.*", "\\1",
    basename(files)
  )))

  for (cid in cond_ids) {
    cid_num <- as.integer(cid)
    fname_pdf <- file.path(
      checks_dir,
      sprintf("checks_cond_%03d.pdf", cid_num)
    )
    pdf(fname_pdf, width = 11, height = 8.5)

    reps <- sort(list.files(data_dir,
      sprintf("sim_data_cond%03d_rep\\d+\\.rds", cid_num),
      full.names = TRUE
    ))
    reps <- sample(reps, min(5, length(reps)))
    for (rf in reps) {
      dat <- readRDS(rf)
      cond_lbl <- sprintf(
        "Cond %03d  (%s)", cid_num,
        dat$true_params$margin1$type
      )
      plot_diagnostics(
        dat$data, dat$true_params,
        cond_lbl, dat$rep_i
      )
    }
    dev.off()
    message("Diagnostics written: ", basename(fname_pdf))
  }
}
