###########################################################################
# simulate_data.R
#
# Generates bivariate VAR(1) time series data with specified marginal
# distributions (Normal or Skew-Normal, potentially different per margin)
# and copulas (Gaussian or Clayton) for the residuals.
# Called by run_pipeline.R
###########################################################################

# --- Libraries ---
if (!requireNamespace("copula", quietly = TRUE)) stop("Package 'copula' is required.")
if (!requireNamespace("sn", quietly = TRUE)) stop("Package 'sn' is required for Skew-Normal.")
if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Package 'mvtnorm' is required.") # Used implicitly by copula pkg

# --- Helper Function: Generate Residuals from Copula ---
generate_residuals <- function(n, margin1_dist, margin1_params,
                               margin2_dist, margin2_params, copula_info) {
  qmargin1 <- switch(margin1_dist,
    "normal" = function(p) qnorm(p, mean = margin1_params$mean, sd = margin1_params$sd),
    "skewnormal" = function(p) sn::qsn(p, xi = margin1_params$xi, omega = margin1_params$omega, alpha = margin1_params$alpha, solver = "RFB"),
    stop("Unsupported margin1 distribution: ", margin1_dist)
  )
  qmargin2 <- switch(margin2_dist,
    "normal" = function(p) qnorm(p, mean = margin2_params$mean, sd = margin2_params$sd),
    "skewnormal" = function(p) sn::qsn(p, xi = margin2_params$xi, omega = margin2_params$omega, alpha = margin2_params$alpha, solver = "RFB"),
    stop("Unsupported margin2 distribution: ", margin2_dist)
  )

  copula_family <- copula_info$type
  copula_param <- copula_info$value

  cop <- tryCatch(
    {
      switch(copula_family,
        "gaussian" = copula::normalCopula(param = copula_param, dim = 2),
        "clayton" = copula::claytonCopula(param = copula_param, dim = 2),
        stop("Unsupported copula type: ", copula_family)
      )
    },
    error = function(e) stop(sprintf("Error creating copula '%s' with param %f: %s", copula_family, copula_param, e$message))
  )

  u <- tryCatch(
    {
      copula::rCopula(n, cop)
    },
    error = function(e) stop(sprintf("Error sampling from copula '%s' with param %f: %s", copula_family, copula_param, e$message))
  )

  residuals <- matrix(NA, nrow = n, ncol = 2)
  colnames(residuals) <- c("e1", "e2")
  eps <- 1e-9
  u[, 1] <- pmax(eps, pmin(1 - eps, u[, 1]))
  u[, 2] <- pmax(eps, pmin(1 - eps, u[, 2]))
  residuals[, 1] <- qmargin1(u[, 1])
  residuals[, 2] <- qmargin2(u[, 2])

  # Check for non-finite values after transformation
  if (any(!is.finite(residuals))) {
    warning(sprintf(
      "Non-finite values generated in residuals (n=%d, m1=%s, m2=%s, cop=%s). Check parameters.",
      n, margin1_dist, margin2_dist, copula_family
    ), call. = FALSE)
    residuals[!is.finite(residuals)] <- NA # Set to NA if problematic
  }
  return(residuals)
}

# --- Function to Simulate One Replication ---
simulate_one_replication_var1 <- function(condition, rep_i, output_dir) {
  cond_id <- condition$condition_id
  T_sim <- condition$T + condition$burnin
  mu <- c(condition$mu1, condition$mu2)
  phi_mat <- matrix(c(condition$phi11, condition$phi21, condition$phi12, condition$phi22), nrow = 2, ncol = 2)

  # Extract margin parameters (handle nesting)
  margin1_params_list <- condition$dgp_margin1_params[[1]]
  margin2_params_list <- condition$dgp_margin2_params[[1]]

  # Generate Residuals
  residuals_all <- tryCatch(
    {
      generate_residuals(
        n = T_sim, margin1_dist = condition$dgp_margin1_dist, margin1_params = margin1_params_list,
        margin2_dist = condition$dgp_margin2_dist, margin2_params = margin2_params_list,
        copula_info = condition$dgp_copula_info[[1]]
      )
    },
    error = function(e) {
      cat(sprintf("\nERROR generating residuals Cond %d, Rep %d: %s\n", cond_id, rep_i, e$message))
      NULL
    }
  )

  if (is.null(residuals_all) || anyNA(residuals_all)) {
    cat(sprintf("  Skipping Rep %d Cond %d: residual generation issues.\n", rep_i, cond_id))
    return(invisible(NULL))
  }

  # Simulate VAR(1) Process
  y <- matrix(NA, nrow = T_sim, ncol = 2)
  colnames(y) <- c("y1", "y2")
  y[1, ] <- c(0, 0) # Start at zero
  for (t in 2:T_sim) {
    y_lagged <- y[t - 1, ]
    if (anyNA(y_lagged)) {
      cat(sprintf("\nERROR: NA in lagged y t=%d Cond %d, Rep %d.\n", t, cond_id, rep_i))
      return(invisible(NULL))
    }
    predicted_mean <- mu + phi_mat %*% y_lagged
    y[t, ] <- predicted_mean + residuals_all[t, ]
  }
  if (any(!is.finite(y))) {
    cat(sprintf("\nWARNING: Non-finite VAR series Cond %d, Rep %d.\n", cond_id, rep_i))
  }

  # Remove Burn-in & Prepare Output
  y_final <- y[(condition$burnin + 1):T_sim, , drop = FALSE]
  T_final <- nrow(y_final)
  sim_data_df <- data.frame(i = 1, t = 1:T_final, y1 = y_final[, 1], y2 = y_final[, 2])

  true_params_list <- list(
    mu = mu, phi = phi_mat,
    margin1_dist = condition$dgp_margin1_dist, margin1_params = margin1_params_list,
    margin2_dist = condition$dgp_margin2_dist, margin2_params = margin2_params_list,
    copula_type = condition$dgp_copula_info[[1]]$type, copula_param_name = condition$dgp_copula_info[[1]]$param_name,
    copula_param_value = condition$dgp_copula_info[[1]]$value, copula_tau = condition$dgp_tau
  )

  output_list <- list(
    condition_id = cond_id, rep_i = rep_i, N = 1, T = T_final, data = sim_data_df,
    true_params = true_params_list,
    dgp_info = list(
      alpha1 = condition$dgp_alpha1, alpha2 = condition$dgp_alpha2,
      copula = condition$dgp_copula_type, tau = condition$dgp_tau
    ),
    correct_fit_model_code = condition$correct_fit_model_code
  )

  # Save Output
  output_filename <- file.path(output_dir, sprintf("sim_data_cond%03d_rep%03d.rds", cond_id, rep_i)) # Simplified name
  tryCatch(
    {
      saveRDS(output_list, file = output_filename)
    },
    error = function(e) {
      cat(sprintf("\nERROR saving file %s: %s\n", output_filename, e$message))
    }
  )

  invisible(NULL)
}

# --- Main Simulation Function ---
simulate_all_conditions_var1 <- function(sim_conditions_df, output_dir) {
  if (!dir.exists(output_dir)) stop("Output directory does not exist: ", output_dir)
  total_conditions <- nrow(sim_conditions_df)
  cat(sprintf("Starting simulation: %d conditions, %d reps each...\n", total_conditions, sim_conditions_df$n_reps[1]))

  for (i in 1:total_conditions) {
    condition <- sim_conditions_df[i, ]
    cond_id <- condition$condition_id
    n_replications <- condition$n_reps
    cat(sprintf(
      "--- Condition ID: %03d (T=%d, A1=%.1f, A2=%.1f, Cop=%s, Tau=%.1f) ---\n",
      cond_id, condition$T, condition$dgp_alpha1, condition$dgp_alpha2, condition$dgp_copula_type, condition$dgp_tau
    ))

    for (j in 1:n_replications) {
      if (j %% max(1, floor(n_replications / 5)) == 0 || j == 1 || j == n_replications) { # Print progress ~5 times
        cat(sprintf("  Running replication %d of %d...\n", j, n_replications))
      }
      simulate_one_replication_var1(condition, j, output_dir)
    }
  }
  cat("\n--- All simulations finished. ---\n")
  invisible(NULL)
}
