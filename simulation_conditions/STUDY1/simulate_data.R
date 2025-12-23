###########################################################################
# simulate_data.R   (single-level, resume-ready version)
#
# Generates one-subject bivariate VAR(1) series with non-Gaussian
# innovations coupled by a Gaussian copula.
# Supports three skewness conditions (moderate / strong / extreme chi-sq).
#
# FIXED: χ² mirror now uses 1-u to preserve copula correlation sign
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {

  stop("Install package 'copula'.")
}
if (!requireNamespace("sn", quietly = TRUE)) {
  stop("Install package 'sn'.")
}

## ========================================================================
## 1. Helper: skew-normal parameters for Var = 1
## ========================================================================
calc_sn_params <- function(alpha, target_var = 1) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(target_var / (1 - 2 * delta^2 / pi))
  xi <- -omega * delta * sqrt(2 / pi)
  list(xi = xi, omega = omega, alpha = alpha)
}

## ========================================================================
## 2. Residual generator (n x 2 matrix)
## ========================================================================
generate_residuals <- function(n, m1, m2, copula_rho, eps = 1e-9) {
  # Helper to draw from copula, resampling boundary values
  draw_u <- function(n, cop) {
    u <- copula::rCopula(n, cop)
    while (any(u <= eps | u >= 1 - eps)) {
      bad <- unique(which(u <= eps | u >= 1 - eps, arr.ind = TRUE)[, 1])
      u[bad, ] <- copula::rCopula(length(bad), cop)
    }
    u
  }

  # Build quantile function for each margin type
  qfun <- function(m) {
    switch(m$type,
      skewnormal = {
        p <- calc_sn_params(m$alpha)
        function(pu) {
          sn::qsn(pu, xi = p$xi, omega = p$omega,
                  alpha = p$alpha, solver = "RFB")
        }
      },
      chisq = {
        # Standardize Chi-Squared (Mean=0, Var=1)
        df <- m$df
        mean_chi <- df
        sd_chi <- sqrt(2 * df)
        mir <- isTRUE(m$mirror)
        
        # CRITICAL FIX: To sample from -X using uniform u while preserving
        # copula correlation, we need Q_{-X}(u) = -Q_X(1-u)
        # Previously we used -Q_X(u) which flips the copula correlation sign!
        function(pu) {
          # Adjust uniform FIRST if mirroring, to preserve copula structure
          pu_adj <- if (mir) 1 - pu else pu
          x_raw <- stats::qchisq(pu_adj, df)
          x_std <- (x_raw - mean_chi) / sd_chi
          # Then negate the result for left-skewed distribution
          if (mir) -x_std else x_std
        }
      },
      stop("Unknown margin type: ", m$type)
    )
  }

  cop <- copula::normalCopula(copula_rho, dim = 2)
  u <- draw_u(n, cop)
  cbind(e1 = qfun(m1)(u[, 1]), e2 = qfun(m2)(u[, 2]))
}

## ========================================================================
## 3. Verify copula correlation (diagnostic helper)
## ========================================================================
verify_copula_correlation <- function(n_test = 10000, m1, m2, copula_rho, tol = 0.05) {
  # Generate test residuals
  eps_test <- generate_residuals(n_test, m1, m2, copula_rho)
  

  # Compute empirical Kendall's tau
  tau_emp <- cor(eps_test[,1], eps_test[,2], method = "kendall")
  
  # Expected tau for Gaussian copula: tau = (2/pi) * arcsin(rho)
 tau_expected <- (2/pi) * asin(copula_rho)
  
  # Check if within tolerance
  tau_diff <- abs(tau_emp - tau_expected)
  
  list(
    tau_empirical = tau_emp,
    tau_expected = tau_expected,
    difference = tau_diff,
    passed = tau_diff < tol
  )
}

## ========================================================================
## 4. Single replication (one design row, one rep)
## ========================================================================
simulate_one_replication <- function(cond_row, rep_i, output_dir, burnin = 100) {
  # Process intercept (0 since innovations are standardized to mean 0)
  mu_vec <- c(0, 0)
  phi_mat <- cond_row$phi_matrix[[1]]
  T_time <- cond_row$T
  T_sim <- T_time + burnin

  eps_mat <- generate_residuals(
    n          = T_sim,
    m1         = cond_row$margin_info[[1]]$margin1,
    m2         = cond_row$margin_info[[1]]$margin2,
    copula_rho = cond_row$rho
  )

  y <- matrix(0, T_sim, 2)
  y[1, ] <- eps_mat[1, ]
  
  for (t in 2:T_sim) {
    y[t, ] <- mu_vec + phi_mat %*% y[t - 1, ] + eps_mat[t, ]
  }

  y_final <- y[(burnin + 1):T_sim, , drop = FALSE]

  out <- list(
    condition_id = cond_row$condition_id,
    rep_i = rep_i,
    N = 1,
    T = T_time,
    data = data.frame(
      i  = 1,
      t  = seq_len(T_time),
      y1 = y_final[, 1],
      y2 = y_final[, 2]
    ),
    phi_matrix = phi_mat,
    rho = cond_row$rho,
    true_params = list(
      mu         = mu_vec,
      phi        = phi_mat,
      copula_rho = cond_row$rho,
      margin1    = cond_row$margin_info[[1]]$margin1,
      margin2    = cond_row$margin_info[[1]]$margin2
    )
  )

  saveRDS(
    out,
    file.path(output_dir,
              sprintf("sim_data_cond%03d_rep%03d.rds",
                      cond_row$condition_id, rep_i))
  )
}

## ========================================================================
## 5. Driver over the design grid (resume-aware)
## ========================================================================
simulate_all_conditions_var1 <- function(sim_conditions_df,
                                         output_dir,
                                         burnin = 100,
                                         start_condition = 1,
                                         start_rep = 1,
                                         overwrite = FALSE,
                                         verify_copula = TRUE) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating single-level data (burn-in = ", burnin, ") ===")

  # Optional: verify copula correlation for each unique margin configuration
 if (verify_copula) {
    message("Verifying copula correlations for each condition type...")
    unique_configs <- unique(sim_conditions_df[, c("skew_level", "direction", "rho")])
    
    for (i in seq_len(nrow(unique_configs))) {
      cfg <- unique_configs[i, ]
      # Find a matching row to get margin_info
      match_row <- sim_conditions_df[
        sim_conditions_df$skew_level == cfg$skew_level &
        sim_conditions_df$direction == cfg$direction &
        sim_conditions_df$rho == cfg$rho, ][1, ]
      
      check <- verify_copula_correlation(
        n_test = 5000,
        m1 = match_row$margin_info[[1]]$margin1,
        m2 = match_row$margin_info[[1]]$margin2,
        copula_rho = match_row$rho
      )
      
      status <- if (check$passed) "OK" else "WARNING"
      message(sprintf("  %s %s rho=%.2f: tau_emp=%.3f, tau_exp=%.3f [%s]",
                      cfg$skew_level, cfg$direction, cfg$rho,
                      check$tau_empirical, check$tau_expected, status))
      
      if (!check$passed) {
        warning(sprintf("Copula correlation mismatch for %s %s rho=%.2f!",
                        cfg$skew_level, cfg$direction, cfg$rho))
      }
    }
  }

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]

    # Resume logic
    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) start_rep else 1

    message(sprintf(
      " Condition %03d  (%s dir=%s T=%d rho=%.2f VAR=%s)",
      cond$condition_id, cond$skew_level, cond$direction,
      cond$T, cond$rho, cond$VARset
    ))

    pb <- utils::txtProgressBar(min = 1, max = cond$n_reps, 
                                 initial = first_rep, style = 3)

    for (r in seq(from = first_rep, to = cond$n_reps)) {
      file_target <- file.path(
        output_dir,
        sprintf("sim_data_cond%03d_rep%03d.rds", cond$condition_id, r)
      )

      if (!overwrite && file.exists(file_target)) {
        utils::setTxtProgressBar(pb, r)
        next
      }

      simulate_one_replication(cond, r, output_dir, burnin)
      utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    cat("\n")
  }
  message("=== Simulation finished ===")
}
