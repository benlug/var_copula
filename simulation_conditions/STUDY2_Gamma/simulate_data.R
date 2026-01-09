###########################################################################
# simulate_data.R
#
# Single-level VAR(1) simulation engine used across studies.
# Supported innovation marginal types:
#   • skewnormal  (Study 1)
#   • chisq       (Study 1)
#   • exponential (Study 2)
#   • gamma       (Study 2-Gamma extension)
#
# All innovations are standardized to mean 0 and variance 1 on each margin.
# A Gaussian copula couples the margins.
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Install package 'copula'.")
}

## ========================================================================
## 1 · helper: skew‑normal parameters (Study 1 compatibility)
## ========================================================================
calc_sn_params <- function(alpha, target_var = 1) {
  if (!requireNamespace("sn", quietly = TRUE)) stop("Package 'sn' required for Skew-Normal.")
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(target_var / (1 - 2 * delta^2 / pi))
  xi <- -omega * delta * sqrt(2 / pi)
  list(xi = xi, omega = omega, alpha = alpha)
}

## ========================================================================
## 2 · residual generator  (n × 2 matrix)
## ========================================================================
#' Generate standardized innovations with specified marginals and Gaussian copula.
#'
#' @param n number of time points
#' @param m1 list describing margin 1 (type + parameters)
#' @param m2 list describing margin 2
#' @param copula_rho Gaussian copula correlation
#' @param eps boundary avoidance for PIT values
#'
#' @return n×2 matrix of standardized innovations
generate_residuals <- function(n, m1, m2, copula_rho, eps = 1e-9) {
  draw_u <- function(n, cop) {
    u <- copula::rCopula(n, cop)
    # Resample any rows that hit the boundary to avoid Inf inv_Phi or q-functions.
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
          sn::qsn(
            pu,
            xi = p$xi, omega = p$omega,
            alpha = p$alpha, solver = "RFB"
          )
        }
      },

      chisq = {
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

      exponential = {
        rate <- m$rate
        mean_exp <- 1 / rate
        sd_exp <- 1 / rate
        mir <- isTRUE(m$mirror)

        function(pu) {
          x_raw <- stats::qexp(pu, rate = rate)
          x_std <- (x_raw - mean_exp) / sd_exp
          if (mir) -x_std else x_std
        }
      },

      gamma = {
        shape <- m$shape
        rate <- m$rate
        if (!is.finite(shape) || shape <= 0) stop("Gamma shape must be > 0")
        if (!is.finite(rate) || rate <= 0) stop("Gamma rate must be > 0")

        mean_gam <- shape / rate
        sd_gam <- sqrt(shape) / rate
        mir <- isTRUE(m$mirror)

        function(pu) {
          x_raw <- stats::qgamma(pu, shape = shape, rate = rate)
          x_std <- (x_raw - mean_gam) / sd_gam
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
## 3 · single replication
## ========================================================================
#' Simulate one VAR(1) replication and save as an RDS.
#'
#' The function is deterministic given seed_sim and the condition row.
simulate_one_replication <- function(cond_row, rep_i,
                                     output_dir, burnin = 30,
                                     seed_sim = NULL) {
  if (!is.null(seed_sim)) set.seed(as.integer(seed_sim))

  mu_vec <- c(0, 0) # standardized innovations
  phi_mat <- cond_row$phi_matrix[[1]]
  T_time <- cond_row$T
  T_sim <- T_time + burnin

  # Effective copula rho after mirroring: mirroring exactly one margin flips sign.
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
      mu              = mu_vec,
      phi             = phi_mat,
      copula_rho       = cond_row$rho,
      copula_rho_input = cond_row$rho,
      copula_rho_eff   = rho_eff,
      margin1         = cond_row$margin_info[[1]]$margin1,
      margin2         = cond_row$margin_info[[1]]$margin2
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

## ========================================================================
## 4 · driver over the design grid
## ========================================================================
#' Simulate all conditions with resume/skip support.
#'
#' @param sim_conditions_df design grid with condition_id and n_reps.
#' @param start_condition first condition_id to simulate
#' @param start_rep first replication within start_condition
#' @param overwrite if TRUE, regenerate existing RDS files
#' @param seed_base base seed used to derive deterministic per-cell seeds
simulate_all_conditions_var1 <- function(sim_conditions_df,
                                         output_dir,
                                         burnin = 30,
                                         start_condition = 1,
                                         start_rep = 1,
                                         overwrite = FALSE,
                                         seed_base = 2026L) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating single‑level data ===")

  # Harmonize DGP label column name for logging.
  if ("dgp_level" %in% names(sim_conditions_df) && !"skew_level" %in% names(sim_conditions_df)) {
    names(sim_conditions_df)[names(sim_conditions_df) == "dgp_level"] <- "skew_level"
  }

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]

    # Resume logic
    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) start_rep else 1

    message(sprintf(
      " Condition %03d  (%s dir=%s T=%d ρ=%.2f VAR=%s)",
      cond$condition_id,
      as.character(cond$skew_level %||% "DGP"),
      as.character(cond$direction),
      as.integer(cond$T),
      as.numeric(cond$rho),
      as.character(cond$VARset)
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

      # Deterministic seed per (condition_id, rep) supports stop/resume.
      seed_sim <- as.integer(seed_base + as.integer(cond$condition_id) * 100000L + as.integer(r))

      simulate_one_replication(cond, r, output_dir, burnin, seed_sim = seed_sim)
      utils::setTxtProgressBar(pb, r)
    }

    close(pb)
    cat("\n")
  }

  message("=== Simulation finished ===")
}
