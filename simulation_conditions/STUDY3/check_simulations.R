###########################################################################
# check_simulations.R
# Visual sanity checks – produces one PDF per condition (max 5 reps each)
# (Compatible with Study 1 and Study 2)
###########################################################################

# Needed base packages only (graphics, grDevices, stats)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Computes residuals based on the TRUE generating parameters (mu and phi)
compute_residuals_TRUE <- function(y, mu, phi) {
  Tn <- nrow(y)
  if (Tn < 2) {
    return(list(res1 = numeric(0), res2 = numeric(0)))
  }
  # VAR(1) definition: y_t = mu + Phi * y_{t-1} + epsilon_t
  # Residuals: epsilon_t = y_t - (mu + Phi * y_{t-1})

  # Create matrix for mu (intercept)
  # Robustness check for mu (assuming 0 if missing/invalid length)
  if (length(mu) != 2) {
    mu <- c(0, 0)
  }
  mu_matrix <- matrix(mu, nrow = Tn - 1, ncol = 2, byrow = TRUE)

  # Calculate predicted values: mu + Phi * y_{t-1}
  # Note: R's matrix multiplication requires t(phi) if y is (T-1)x2
  predicted <- mu_matrix + y[1:(Tn - 1), ] %*% t(phi)

  res <- y[2:Tn, ] - predicted
  list(res1 = res[, 1], res2 = res[, 2])
}

plot_diagnostics <- function(df, true_par, cond_lbl, rep_lbl) {
  ymat <- as.matrix(df[, c("y1", "y2")])
  # Calculate residuals using the true parameters stored during simulation
  res <- compute_residuals_TRUE(ymat, true_par$mu, true_par$phi)

  # Set up plotting layout
  old <- par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old))

  # Time series plots
  plot(df$t, df$y1, type = "l", main = "Series y1", xlab = "t", ylab = "y1")
  plot(df$t, df$y2, type = "l", main = "Series y2", xlab = "t", ylab = "y2")
  plot(df$y1, df$y2,
    pch = 19, col = "#00000044",
    main = "Scatter y1 vs y2", xlab = "y1", ylab = "y2"
  )

  # Residual diagnostics
  if (length(res$res1) > 0) {
    # Use 'fd' for robust break calculation
    breaks1 <- tryCatch(hist(res$res1, breaks = "fd", plot = FALSE)$breaks, error = function(e) 20)
    breaks2 <- tryCatch(hist(res$res2, breaks = "fd", plot = FALSE)$breaks, error = function(e) 20)

    # Add skew direction info to histogram titles if available (Study 2: Exponential)
    m1_mirror <- true_par$margin1$mirror %||% NA
    m2_mirror <- true_par$margin2$mirror %||% NA

    # Determine direction label: (Left) if mirrored, (Right) otherwise, empty if NA
    get_dir_label <- function(mirror_flag) {
      if (is.na(mirror_flag)) {
        return("")
      }
      if (mirror_flag) {
        return("(Left)")
      } else {
        return("(Right)")
      }
    }

    title1 <- paste("Hist Resid y1", get_dir_label(m1_mirror))
    title2 <- paste("Hist Resid y2", get_dir_label(m2_mirror))

    hist(res$res1, breaks = breaks1, main = title1, xlab = "res1")
    hist(res$res2, breaks = breaks2, main = title2, xlab = "res2")

    # QQ-Norm plots to check deviation from normality
    qqnorm(res$res1, main = "QQ-Norm y1")
    qqline(res$res1)
    qqnorm(res$res2, main = "QQ-Norm y2")
    qqline(res$res2)

    # Autocorrelation of residuals (should be white noise)
    acf(res$res1, main = "ACF res1", na.action = na.pass)
    acf(res$res2, main = "ACF res2", na.action = na.pass)
  }


  # Main Title
  mtext(sprintf("%s  |  Rep %d", cond_lbl, rep_lbl),
    outer = TRUE, line = 0, cex = 1.1, font = 2
  )
}

run_post_sim_checks_var1 <- function(data_dir, checks_dir) {
  files <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(files)) {
    message("No simulation files found in ", data_dir)
    return(invisible())
  }

  cond_ids <- sort(unique(sub(
    "^.*_cond(\\d+)_rep.*", "\\1",
    basename(files)
  )))

  message("Generating diagnostic plots for ", length(cond_ids), " conditions in ", checks_dir)
  dir.create(checks_dir, FALSE, TRUE) # Ensure directory exists

  for (cid in cond_ids) {
    cid_num <- as.integer(cid)
    fname_pdf <- file.path(
      checks_dir,
      sprintf("checks_cond_%03d.pdf", cid_num)
    )


    # Find all replications for this condition
    reps <- list.files(data_dir,
      sprintf("sim_data_cond%03d_rep\\d+\\.rds", cid_num),
      full.names = TRUE
    )

    # Sample up to 5 replications to check
    if (length(reps) > 5) {
      # Set seed for reproducible sampling
      set.seed(cid_num)
      reps_to_check <- sample(reps, 5)
    } else {
      reps_to_check <- reps
    }

    if (length(reps_to_check) > 0) {
      # Use tryCatch around PDF generation for robustness
      tryCatch({
        pdf(fname_pdf, width = 11, height = 8.5)

        for (rf in reps_to_check) {
          dat <- try(readRDS(rf), silent = TRUE)
          if (inherits(dat, "try-error")) {
            message("Error reading RDS file: ", rf)
            next
          }

          # Create condition label - dynamically extracting the distribution type
          dist_type <- dat$true_params$margin1$type %||% "Unknown"

          cond_lbl <- sprintf(
            "Cond %03d (%s)", cid_num, dist_type
          )

          plot_diagnostics(
            dat$data, dat$true_params,
            cond_lbl, dat$rep_i
          )
        }
      }, error = function(e) {
        message("Error generating PDF ", fname_pdf, ": ", e$message)
      }, finally = {
        # Ensure the graphics device is closed even if an error occurs
        if (names(dev.cur()) != "null device") dev.off()
      })

      # message("Diagnostics written: ", basename(fname_pdf))
    }
  }
  message("Simulation checks complete.")
}
