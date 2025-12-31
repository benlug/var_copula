###########################################################################
# simulate_data.R   (Updated for Study II: Exponential DGP)
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Install package 'copula'.")
}
# 'sn' is required if Study 1 (Skew-Normal) functionality is retained.
# if (!requireNamespace("sn", quietly = TRUE)) {
#   stop("Install package 'sn'.")
# }

## ========================================================================
## 1 · helper: skew‑normal parameters (Keep for Study 1 compatibility)
## ========================================================================
calc_sn_params <- function(alpha, target_var = 1) {
  if (!requireNamespace("sn", quietly = TRUE)) stop("Package 'sn' required for Skew-Normal.")
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
        # (Study 1 implementation)
        df <- m$df
        mean_chi <- df
        sd_chi <- sqrt(2 * df)
        mir <- isTRUE(m$mirror)

        function(pu) {
          x_raw <- stats::qchisq(pu, df)
          x_std <- (x_raw - mean_chi) / sd_chi
          if (mir) -x_std else x_std
        }
      },
      # STUDY 2 ADDITION: Standardized Exponential
      exponential = {
        rate <- m$rate # Expected to be 1 in the design
        # Mean and SD of Exponential(lambda) are both 1/lambda
        mean_exp <- 1 / rate
        sd_exp <- 1 / rate
        mir <- isTRUE(m$mirror)

        function(pu) {
          # 1. Generate raw exponential quantile
          x_raw <- stats::qexp(pu, rate)
          # 2. Standardize (Shift to Mean 0, Scale to SD 1)
          x_std <- (x_raw - mean_exp) / sd_exp
          # 3. Mirror if needed (for left skew)
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
## 3 · single replication (No changes needed)
## ========================================================================
simulate_one_replication <- function(cond_row, rep_i,
                                     output_dir, burnin = 30,
                                     seed_sim = NULL) {
  # Deterministic per-dataset seeding (critical for reproducibility under resume/skip)
  if (!is.null(seed_sim)) {
    set.seed(as.integer(seed_sim))
  }
  # mu_vec is 0 as innovations are standardized.
  mu_vec <- c(0, 0)
  phi_mat <- cond_row$phi_matrix[[1]]
  T_time <- cond_row$T
  T_sim <- T_time + burnin

  # Effective copula rho on the PIT/latent-normal scale after mirroring.
  # Mirroring exactly one margin flips the sign of the Gaussian copula correlation.
  mirror1 <- isTRUE(cond_row$margin_info[[1]]$margin1$mirror)
  mirror2 <- isTRUE(cond_row$margin_info[[1]]$margin2$mirror)
  s1 <- if (mirror1) -1 else 1
  s2 <- if (mirror2) -1 else 1
  rho_eff <- s1 * s2 * cond_row$rho

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
    seed_sim = seed_sim,
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
    rho_eff = rho_eff,
    true_params = list(
      mu         = mu_vec,
      phi        = phi_mat,
      # Backward-compatible names (some downstream scripts may still expect `copula_rho`)
      copula_rho       = cond_row$rho,
      copula_rho_input = cond_row$rho,
      copula_rho_eff   = rho_eff,
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
## 4 · driver over the design grid (Minor update for logging)
## ========================================================================
simulate_all_conditions_var1 <- function(sim_conditions_df,
                                         output_dir,
                                         burnin = 30,
                                         start_condition = 1,
                                         start_rep = 1,
                                         overwrite = FALSE,
                                         seed_base = 2026L) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating single‑level data ===")

  # Handle different column names between Study 1 (skew_level) and Study 2 (dgp_level)
  # This makes the logging consistent.
  if ("dgp_level" %in% names(sim_conditions_df) && !"skew_level" %in% names(sim_conditions_df)) {
    names(sim_conditions_df)[names(sim_conditions_df) == "dgp_level"] <- "skew_level"
  }

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

      # Deterministic seed per (condition, replication) ensures reproducibility under resume/skip
      # Use a large multiplier to avoid collisions if reps/conditions are expanded later.
      seed_sim <- as.integer(seed_base + as.integer(cond$condition_id) * 100000L + as.integer(r))

      simulate_one_replication(cond, r, output_dir, burnin, seed_sim = seed_sim)
      utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    cat("\n")
  }
  message("=== Simulation finished ===")
}
