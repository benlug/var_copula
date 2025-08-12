###########################################################################
# simulate_data.R   (single‑level, resume‑ready version)
#
# · Generates one‑subject bivariate VAR(1) series with non‑Gaussian
#   innovations coupled by a Gaussian copula.
# · Supports three skewness conditions (moderate / strong / extreme χ²¹).
# · Called from run_pipeline.R, which passes the design grid.
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Install package 'copula'.")
}
if (!requireNamespace("sn", quietly = TRUE)) {
  stop("Install package 'sn'.")
}

## ========================================================================
## 1 · helper: skew‑normal parameters, Var = 1
## ========================================================================
calc_sn_params <- function(alpha, target_var = 1) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(target_var / (1 - 2 * delta^2 / pi))
  xi <- -omega * delta * sqrt(2 / pi)
  list(xi = xi, omega = omega, alpha = alpha)
}

## ========================================================================
## 2 · residual generator  (n × 2 matrix)
## ========================================================================
generate_residuals <- function(n, m1, m2, copula_rho, eps = 1e-9) {
  draw_u <- function(n, cop) { # resample p ≈ 0/1 rows
    u <- copula::rCopula(n, cop)
    while (any(u <= eps | u >= 1 - eps)) {
      bad <- unique(which(u <= eps | u >= 1 - eps, arr.ind = TRUE)[, 1])
      u[bad, ] <- copula::rCopula(length(bad), cop)
    }
    u
  }

  qfun <- function(m) {
    switch(m$type,
      skewnormal = {
        p <- calc_sn_params(m$alpha)
        function(pu) {
          sn::qsn(pu,
            xi = p$xi, omega = p$omega,
            alpha = p$alpha, solver = "RFB"
          )
        }
      },
      chisq = {
        # FIX 1: Standardize Chi-Squared (Mean=0, Var=1)
        df <- m$df
        mean_chi <- df
        sd_chi <- sqrt(2 * df)
        mir <- isTRUE(m$mirror)
        
        function(pu) {
          x_raw <- stats::qchisq(pu, df)
          # Standardize
          x_std <- (x_raw - mean_chi) / sd_chi
          # Mirror if needed (standardization ensures mirroring around 0)
          if (mir) -x_std else x_std
        }
      },
      stop("unknown margin type")
    )
  }

  cop <- copula::normalCopula(copula_rho, dim = 2)
  u <- draw_u(n, cop)
  cbind(e1 = qfun(m1)(u[, 1]), e2 = qfun(m2)(u[, 2]))
}

## ========================================================================
## 3 · single replication (one design row · one rep)
## ========================================================================
simulate_one_replication <- function(cond_row, rep_i,
                                     output_dir, burnin = 30) {
  # Note: mu_vec is the process intercept. Since innovations are standardized 
  # (Mean=0, Var=1) for all conditions now, the process mean is also 0.
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
  # Initialize first time point using the first residual for better stability
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
    file.path(
      output_dir,
      sprintf(
        "sim_data_cond%03d_rep%03d.rds",
        cond_row$condition_id, rep_i
      )
    )
  )
}

## ========================================================================
## 4 · driver over the design grid  (resume‑aware)
## ========================================================================
simulate_all_conditions_var1 <- function(sim_conditions_df,
                                         output_dir,
                                         burnin = 30,
                                         start_condition = 1,
                                         start_rep = 1,
                                         overwrite = FALSE) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating single‑level data ===")

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]

    ## -- resume logic ----------------------------------------------------
    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) {
      start_rep
    } else {
      1
    }

    message(sprintf(
      " Condition %03d  (%s dir=%s T=%d ρ=%.2f VAR=%s)",
      cond$condition_id, cond$skew_level, cond$direction,
      cond$T, cond$rho, cond$VARset
    ))

    # FIX 2: Correct progress bar initialization
    # Initialize min=1 (or 0) and set the starting point with 'initial'
    pb <- utils::txtProgressBar(min = 1, max = cond$n_reps, initial = first_rep, style = 3)

    for (r in seq(from = first_rep, to = cond$n_reps)) {
      file_target <- file.path(
        output_dir,
        sprintf(
          "sim_data_cond%03d_rep%03d.rds",
          cond$condition_id, r
        )
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