###########################################################################
# simulate_data_tvp.R   (Study 6: Time-Varying Parameter Copula-VAR)
###########################################################################
#
# Generates bivariate VAR(1) time series with:
#   - Optional time-varying copula correlation rho_t
#   - Various rho trajectory patterns (constant, linear, step, oscillating, U-shaped)
#   - Gaussian or Exponential marginal innovations
#
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Install package 'copula'.")
}

## ========================================================================
## 1. Generate rho trajectory based on pattern type
## ========================================================================

generate_rho_trajectory <- function(T, pattern, params) {
  # Returns a vector of length T with rho values for each time point
  #
  # pattern: "constant", "linear", "step", "oscillating", "u_shaped", "random_walk"
  # params: list with pattern-specific parameters
  
  t_seq <- 1:T
  t_norm <- t_seq / T  # Normalized time [0, 1]
  
  rho_t <- switch(pattern,
    
    # Pattern 1: Constant (baseline/null case)
    constant = {
      rep(params$rho_base, T)
    },
    
    # Pattern 2: Linear drift (e.g., therapy reducing coupling)
    linear = {
      rho_start <- params$rho_start
      rho_end <- params$rho_end
      rho_start + (rho_end - rho_start) * t_norm
    },
    
    # Pattern 3: Step change (e.g., episode onset)
    step = {
      rho_before <- params$rho_before
      rho_after <- params$rho_after
      tau <- params$tau  # Change point (proportion of T)
      ifelse(t_norm <= tau, rho_before, rho_after)
    },
    
    # Pattern 4: Oscillating (e.g., weekly cycles)
    oscillating = {
      rho_mean <- params$rho_mean
      amplitude <- params$amplitude
      period <- params$period  # In time points
      rho_mean + amplitude * sin(2 * pi * t_seq / period)
    },
    
    # Pattern 5: U-shaped (temporary spike during crisis)
    u_shaped = {
      rho_baseline <- params$rho_baseline
      delta_rho <- params$delta_rho
      # Quadratic: peaks at t_norm = 0.5
      rho_baseline + delta_rho * (1 - 4 * (t_norm - 0.5)^2)
    },
    
    # Pattern 6: Random walk (most general TVP)
    random_walk = {
      rho_start <- params$rho_start
      sigma_z <- params$sigma_z
      
      # Generate on z-scale (Fisher transformation), then transform
      z <- numeric(T)
      z[1] <- atanh(rho_start)  # Initial state
      
      for (t in 2:T) {
        z[t] <- z[t-1] + rnorm(1, 0, sigma_z)
      }
      
      # Transform back to rho scale
      tanh(z)
    },
    
    stop("Unknown pattern type: ", pattern)
  )
  
  # Ensure rho stays within valid bounds
  rho_t <- pmax(-0.99, pmin(0.99, rho_t))
  
  return(rho_t)
}


## ========================================================================
## 2. Residual generator with time-varying copula
## ========================================================================

generate_residuals_tvp <- function(T, m1, m2, rho_trajectory, eps = 1e-9) {
  # Generate residuals with time-varying copula correlation
  #
  # T: Number of time points
  # m1, m2: Marginal distribution specifications
  # rho_trajectory: Vector of length T with rho_t values
  
  if (length(rho_trajectory) != T) {
    stop("rho_trajectory must have length T")
  }
  
  # Quantile function for each margin
  qfun <- function(m) {
    switch(m$type,
      normal = {
        function(pu) qnorm(pu, 0, 1)
      },
      exponential = {
        rate <- m$rate
        mean_exp <- 1 / rate
        sd_exp <- 1 / rate
        mir <- isTRUE(m$mirror)
        
        function(pu) {
          x_raw <- qexp(pu, rate)
          x_std <- (x_raw - mean_exp) / sd_exp
          if (mir) -x_std else x_std
        }
      },
      stop("Unknown margin type: ", m$type)
    )
  }
  
  q1 <- qfun(m1)
  q2 <- qfun(m2)
  
  # Generate residuals for each time point with time-specific rho
  e1 <- numeric(T)
  e2 <- numeric(T)
  
  for (t in 1:T) {
    rho_t <- rho_trajectory[t]
    
    # Generate from bivariate normal with correlation rho_t
    # Using Cholesky decomposition
    w1 <- rnorm(1)
    w2 <- rho_t * w1 + sqrt(1 - rho_t^2) * rnorm(1)
    
    # Transform to uniform via normal CDF
    u1 <- pnorm(w1)
    u2 <- pnorm(w2)
    
    # Clamp to avoid boundary issues
    u1 <- max(eps, min(1 - eps, u1))
    u2 <- max(eps, min(1 - eps, u2))
    
    # Transform to marginal scale
    e1[t] <- q1(u1)
    e2[t] <- q2(u2)
  }
  
  cbind(e1 = e1, e2 = e2)
}


## ========================================================================
## 3. Simulate one replication with TVP
## ========================================================================

simulate_one_replication_tvp <- function(cond_row, rep_i, output_dir, burnin = 50) {
  
  mu_vec <- c(0, 0)
  phi_mat <- cond_row$phi_matrix[[1]]
  T_time <- cond_row$T
  T_sim <- T_time + burnin
  
  # Generate rho trajectory for simulation period
  rho_traj_full <- generate_rho_trajectory(
    T = T_sim,
    pattern = cond_row$tvp_pattern,
    params = cond_row$tvp_params[[1]]
  )
  
  # Generate residuals with time-varying copula
  eps_mat <- generate_residuals_tvp(
    T = T_sim,
    m1 = cond_row$margin_info[[1]]$margin1,
    m2 = cond_row$margin_info[[1]]$margin2,
    rho_trajectory = rho_traj_full
  )
  
  # Simulate VAR process
  y <- matrix(0, T_sim, 2)
  y[1, ] <- eps_mat[1, ]
  
  for (t in 2:T_sim) {
    y[t, ] <- mu_vec + phi_mat %*% y[t - 1, ] + eps_mat[t, ]
  }
  
  # Discard burn-in
  y_final <- y[(burnin + 1):T_sim, , drop = FALSE]
  rho_traj_final <- rho_traj_full[(burnin + 1):T_sim]
  
  # Compute effective rho for mirrored margins
  m1_mirror <- cond_row$margin_info[[1]]$margin1$mirror %||% FALSE
  m2_mirror <- cond_row$margin_info[[1]]$margin2$mirror %||% FALSE
  s1 <- ifelse(m1_mirror, -1, 1)
  s2 <- ifelse(m2_mirror, -1, 1)
  rho_eff_traj <- s1 * s2 * rho_traj_final
  
  # Output structure
  out <- list(
    condition_id = cond_row$condition_id,
    rep_i = rep_i,
    N = 1,
    T = T_time,
    data = data.frame(
      i = 1,
      t = seq_len(T_time),
      y1 = y_final[, 1],
      y2 = y_final[, 2]
    ),
    phi_matrix = phi_mat,
    true_params = list(
      mu = mu_vec,
      phi = phi_mat,
      rho_trajectory = rho_traj_final,
      rho_eff_trajectory = rho_eff_traj,
      rho_mean = mean(rho_traj_final),
      rho_sd = sd(rho_traj_final),
      rho_start = rho_traj_final[1],
      rho_end = rho_traj_final[T_time],
      rho_range = max(rho_traj_final) - min(rho_traj_final),
      tvp_pattern = cond_row$tvp_pattern,
      tvp_params = cond_row$tvp_params[[1]],
      margin1 = cond_row$margin_info[[1]]$margin1,
      margin2 = cond_row$margin_info[[1]]$margin2
    )
  )
  
  saveRDS(
    out,
    file.path(
      output_dir,
      sprintf("sim_data_cond%03d_rep%03d.rds", cond_row$condition_id, rep_i)
    )
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ========================================================================
## 4. Driver over design grid
## ========================================================================

simulate_all_conditions_tvp <- function(sim_conditions_df,
                                        output_dir,
                                        burnin = 50,
                                        start_condition = 1,
                                        start_rep = 1,
                                        overwrite = FALSE) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating TVP Copula-VAR data ===")
  
  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]
    
    # Resume logic
    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) start_rep else 1
    
    message(sprintf(
      " Condition %03d  (Pattern=%s T=%d margin=%s direction=%s VARset=%s)",
      cond$condition_id, cond$tvp_pattern, cond$T,
      cond$margin_type, cond$direction, cond$VARset
    ))
    
    pb <- utils::txtProgressBar(min = 1, max = cond$n_reps, initial = first_rep, style = 3)
    
    for (r in seq(from = first_rep, to = cond$n_reps)) {
      file_target <- file.path(
        output_dir,
        sprintf("sim_data_cond%03d_rep%03d.rds", cond$condition_id, r)
      )
      
      if (!overwrite && file.exists(file_target)) {
        utils::setTxtProgressBar(pb, r)
        next
      }
      
      simulate_one_replication_tvp(cond, r, output_dir, burnin)
      utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    cat("\n")
  }
  message("=== TVP Simulation finished ===")
}
