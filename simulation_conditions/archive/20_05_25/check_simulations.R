###########################################################################
# check_simulations.R
#
# Performs visual diagnostic checks on simulated single-level VAR(1) data.
# Generates multi-panel plots for each replication/condition.
# Creates one PDF per condition containing plots.
# Called by run_pipeline.R
###########################################################################

# --- Libraries ---
# Assumes dplyr, grDevices, stats, graphics are available

# --- Helper Functions ---

compute_residuals_trueVAR <- function(y, mu, phi_mat) {
  # compute residuals under the true var(1) model
  if (!is.numeric(mu) || length(mu) != 2) stop("mu must be numeric vector length 2.")
  if (!is.matrix(phi_mat) || !all(dim(phi_mat) == c(2, 2))) stop("phi_mat must be 2x2 matrix.")
  if (!is.matrix(y)) y <- as.matrix(y)
  if (ncol(y) != 2) stop("Input 'y' must have 2 columns.")
  T_val <- nrow(y)
  if (T_val < 2) {
    return(list(res1 = numeric(0), res2 = numeric(0)))
  }

  res1 <- numeric(T_val - 1)
  res2 <- numeric(T_val - 1)
  for (t in 2:T_val) {
    y_lagged <- y[t - 1, ]
    # var(1) mean equation
    cond_mean <- mu + phi_mat %*% y_lagged
    residuals_t <- y[t, ] - cond_mean
    res1[t - 1] <- residuals_t[1]
    res2[t - 1] <- residuals_t[2]
  }
  list(res1 = res1, res2 = res2)
}

plot_checks_for_replication <- function(df_rep, true_params, dgp_info, cond_id, rep_id) {
  # produce diagnostic plots for one replication
  acf_fun <- stats::acf
  pacf_fun <- stats::pacf
  if (requireNamespace("forecast", quietly = TRUE)) {
    acf_fun <- forecast::Acf
    pacf_fun <- forecast::Pacf
  }
  plot_error_message <- function(msg) {
    plot.new()
    tryCatch(title(main = msg, col.main = "red", cex.main = 0.8), error = function(e) {})
    box()
  }

  if (!is.data.frame(df_rep) || nrow(df_rep) < 2 || !all(c("t", "y1", "y2") %in% names(df_rep))) {
    plot_error_message(sprintf("Rep %d\nInvalid Data", rep_id))
    return()
  }
  y_mat <- as.matrix(df_rep[, c("y1", "y2")])

  if (is.null(true_params$mu) || is.null(true_params$phi)) {
    plot_error_message(sprintf("Rep %d\nNo True Params", rep_id))
    return()
  }
  res_list <- tryCatch(compute_residuals_trueVAR(y_mat, true_params$mu, true_params$phi), error = function(e) NULL)
  if (is.null(res_list) || length(res_list$res1) == 0) {
    plot_error_message(sprintf("Rep %d\nResid Calc Fail", rep_id))
    return()
  }
  res1 <- res_list$res1
  res2 <- res_list$res2

  old_par <- par(mfrow = c(4, 4), mar = c(3, 3.5, 2.5, 1), oma = c(0, 0, 3, 0), mgp = c(2, 0.7, 0))
  on.exit(par(old_par), add = TRUE)

  # Row 1: Time Series
  tryCatch(plot(df_rep$t, y_mat[, 1], type = "l", col = "blue", main = "Time Series y1", xlab = "Time", ylab = "y1"), error = function(e) plot_error_message("Plot y1 Err"))
  tryCatch(plot(df_rep$t, y_mat[, 2], type = "l", col = "red", main = "Time Series y2", xlab = "Time", ylab = "y2"), error = function(e) plot_error_message("Plot y2 Err"))
  plot.new()
  box()
  plot.new()
  box() # Spacers

  # Row 2: Scatter & Residual Histograms
  tryCatch(plot(y_mat[, 1], y_mat[, 2], main = "Scatter (y1 vs y2)", xlab = "y1", ylab = "y2", pch = 19, cex = 0.7, col = adjustcolor("purple", alpha.f = 0.5)), error = function(e) plot_error_message("Scatter Err"))
  tryCatch(hist(res1, breaks = "Sturges", main = "Hist True Resid (y1)", xlab = "Residual", col = "lightblue"), error = function(e) plot_error_message("Hist R1 Err"))
  tryCatch(hist(res2, breaks = "Sturges", main = "Hist True Resid (y2)", xlab = "Residual", col = "lightcoral"), error = function(e) plot_error_message("Hist R2 Err"))
  plot.new()
  box()

  # Row 3: Residual QQ Plots & Time Plots
  tryCatch(
    {
      qqnorm(res1, main = "QQ True Resid (y1)", pch = 19, cex = 0.6)
      qqline(res1, col = "red", lwd = 1.5)
    },
    error = function(e) plot_error_message("QQ R1 Err")
  )
  tryCatch(
    {
      qqnorm(res2, main = "QQ True Resid (y2)", pch = 19, cex = 0.6)
      qqline(res2, col = "blue", lwd = 1.5)
    },
    error = function(e) plot_error_message("QQ R2 Err")
  )
  tryCatch(plot(seq_along(res1), res1, type = "l", main = "Time Plot True Resid (y1)", xlab = "Index", ylab = "Resid", col = "grey40"), error = function(e) plot_error_message("Plot R1 Err"))
  abline(h = 0, col = "red", lty = 2)
  tryCatch(plot(seq_along(res2), res2, type = "l", main = "Time Plot True Resid (y2)", xlab = "Index", ylab = "Resid", col = "grey40"), error = function(e) plot_error_message("Plot R2 Err"))
  abline(h = 0, col = "blue", lty = 2)

  # Row 4: ACF/PACF of Original Series
  tryCatch(acf_fun(y_mat[, 1], lag.max = 20, main = "ACF of y1", na.action = na.pass), error = function(e) plot_error_message("ACF y1 Err"))
  tryCatch(pacf_fun(y_mat[, 1], lag.max = 20, main = "PACF of y1", na.action = na.pass), error = function(e) plot_error_message("PACF y1 Err"))
  tryCatch(acf_fun(y_mat[, 2], lag.max = 20, main = "ACF of y2", na.action = na.pass), error = function(e) plot_error_message("ACF y2 Err"))
  tryCatch(pacf_fun(y_mat[, 2], lag.max = 20, main = "PACF of y2", na.action = na.pass), error = function(e) plot_error_message("PACF y2 Err"))

  # Updated Title
  title_text <- sprintf(
    "Cond %s Rep %s | DGP: A1=%.1f, A2=%.1f, Cop=%s(tau=%.1f)",
    cond_id, rep_id, dgp_info$alpha1, dgp_info$alpha2, dgp_info$copula, dgp_info$tau
  )
  tryCatch(mtext(title_text, outer = TRUE, line = 1, cex = 1.1, font = 2), error = function(e) {})
  invisible(NULL)
}

# --- Main Function ---
run_post_sim_checks_var1 <- function(data_dir, checks_dir) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' needed.")
  if (!requireNamespace("grDevices", quietly = TRUE)) stop("Package 'grDevices' needed.")
  if (!dir.exists(data_dir)) stop("Data directory not found: ", data_dir)
  if (!dir.exists(checks_dir)) dir.create(checks_dir, recursive = TRUE)

  # Use simplified data file name pattern
  data_files <- list.files(data_dir, pattern = "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (length(data_files) == 0) {
    cat(sprintf("No simulation data files found in '%s'.\n", data_dir))
    return(invisible(NULL))
  }

  cond_ids_chr <- unique(sub("^.*_cond(\\d+)_rep\\d+\\.rds$", "\\1", basename(data_files)))
  cond_ids <- sort(as.integer(cond_ids_chr[!is.na(as.integer(cond_ids_chr))]))
  if (length(cond_ids) == 0) {
    cat(sprintf("Could not extract condition IDs from '%s'.\n", data_dir))
    return(invisible(NULL))
  }
  cat(sprintf("Found data for %d conditions.\n", length(cond_ids)))

  for (cond_id in cond_ids) {
    cat(sprintf("Processing Condition ID %d...\n", cond_id))
    cond_pattern <- sprintf("_cond%03d_rep", cond_id)
    cond_files <- data_files[grepl(cond_pattern, basename(data_files))]
    if (length(cond_files) == 0) {
      cat(sprintf("  No files found for Cond %d.\n", cond_id))
      next
    }

    # Simplified PDF filename
    # save plots for this condition
    pdf_file <- file.path(checks_dir, sprintf("checks_cond_%03d.pdf", cond_id))
    pdf_opened <- FALSE
    tryCatch(
      {
        grDevices::pdf(pdf_file, width = 16, height = 10)
        pdf_opened <- TRUE
      },
      error = function(e) cat(sprintf("  ERROR opening PDF %s: %s\n", pdf_file, e$message))
    )
    if (!pdf_opened) next

    reps_to_plot <- sample(cond_files, min(length(cond_files), 5)) # Plot max 5 reps
    cat(sprintf("  Plotting %d replication(s) for condition %d...\n", length(reps_to_plot), cond_id))

    for (f in reps_to_plot) {
      rep_pattern <- sprintf("(?<=_cond%03d_rep)(\\d+)(?=\\.rds$)", cond_id)
      rep_id_match <- regmatches(basename(f), regexpr(rep_pattern, basename(f), perl = TRUE))
      if (length(rep_id_match) == 0) {
        cat("    Cannot extract rep number from:", basename(f), "\n")
        next
      }
      rep_id <- as.integer(rep_id_match)
      # cat(sprintf("    Plotting Rep %d...\n", rep_id)) # Reduced verbosity

      sim_dat <- tryCatch(
        {
          readRDS(f)
        },
        error = function(e) {
          cat("      ERROR reading:", e$message, "\n")
          NULL
        }
      )
      if (is.null(sim_dat)) next

      # Updated Validation Check for alpha1/alpha2
      valid_structure <- is.list(sim_dat) && !is.null(sim_dat$data) && is.data.frame(sim_dat$data) &&
        !is.null(sim_dat$true_params) && is.list(sim_dat$true_params) &&
        !is.null(sim_dat$dgp_info) && is.list(sim_dat$dgp_info) &&
        all(c("t", "y1", "y2") %in% names(sim_dat$data)) &&
        !is.null(sim_dat$true_params$mu) && length(sim_dat$true_params$mu) == 2 &&
        !is.null(sim_dat$true_params$phi) && is.matrix(sim_dat$true_params$phi) && all(dim(sim_dat$true_params$phi) == c(2, 2)) &&
        !is.null(sim_dat$dgp_info$alpha1) && !is.null(sim_dat$dgp_info$alpha2) && # Check for both alphas
        !is.null(sim_dat$dgp_info$copula) && !is.null(sim_dat$dgp_info$tau)

      if (!valid_structure) {
        cat("      WARNING: Skipping rep", rep_id, "- invalid structure.\n")
        next
      }

      tryCatch(
        {
          plot_checks_for_replication(
            df_rep = sim_dat$data, true_params = sim_dat$true_params,
            dgp_info = sim_dat$dgp_info, cond_id = cond_id, rep_id = rep_id
          )
        },
        error = function(e) {
          cat(sprintf("      ERROR plotting Rep %d: %s\n", rep_id, e$message))
        }
      )
    } # End replication loop

    tryCatch(grDevices::dev.off(), error = function(e) { }) # Silence error if already closed
    cat(sprintf("  Finished plots => %s\n", pdf_file))

  } # End condition loop

  cat(sprintf("\n--- Simulation checks finished. Output PDFs in '%s'. ---\n", checks_dir))
  invisible(NULL)
}
