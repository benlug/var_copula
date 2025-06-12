###########################################################################
# simulate_data.R   (single‑level, 192‑cell design)
#
# • Supports Skew‑Normal (±4, ±9) and χ²(df=1) residuals
# • χ² left‑skew created via sign‑mirroring
# • Uses Gaussian copula with target ρ supplied by design grid
# • Autoregressive parameters read from design row
###########################################################################

# Required CRAN packages:
#   copula     (for normalCopula)
#   sn         (for skew‑normal quantiles)
#   mvtnorm    (via copula)

if (!requireNamespace("copula", quietly = TRUE)) stop("Install package 'copula'.")
if (!requireNamespace("sn", quietly = TRUE)) stop("Install package 'sn'.")

# ---------- helper to choose skew‑normal parameters with Var=1 ------------
calc_sn_params <- function(alpha, target_var = 1) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(target_var / (1 - (2 * delta^2) / pi))
  xi <- -omega * delta * sqrt(2 / pi)
  list(xi = xi, omega = omega, alpha = alpha)
}

# ------------------------------------------------------------------------
#  Robust residual generator (Skew‑Normal | χ², Gaussian copula)
# ------------------------------------------------------------------------
generate_residuals <- function(n, m1, m2, copula_rho, eps = 1e-9) {
  # --- helper: redraw rows that fall too close to 0 or 1 -----------------
  draw_u <- function(n, cop, eps) {
    u <- copula::rCopula(n, cop) # initial draw
    # Any rows with problematic values?
    bad <- which(u <= eps | u >= 1 - eps, arr.ind = TRUE)
    if (nrow(bad) == 0L) {
      return(u)
    }
    repeat {
      rows_bad <- unique(bad[, 1])
      u[rows_bad, ] <- copula::rCopula(length(rows_bad), cop)
      bad <- which(u <= eps | u >= 1 - eps, arr.ind = TRUE)
      if (nrow(bad) == 0L) break
    }
    u
  }

  # --- build marginal quantile functions --------------------------------
  q1 <- switch(m1$type,
    skewnormal = {
      par <- calc_sn_params(m1$alpha)
      function(p) {
        sn::qsn(p,
          xi = par$xi, omega = par$omega,
          alpha = par$alpha, solver = "RFB"
        )
      }
    },
    chisq = {
      df <- m1$df
      mir <- isTRUE(m1$mirror)
      function(p) {
        x <- stats::qchisq(p, df = df)
        if (mir) -x else x
      }
    },
    stop("unknown margin1 type")
  )

  q2 <- switch(m2$type,
    skewnormal = {
      par <- calc_sn_params(m2$alpha)
      function(p) {
        sn::qsn(p,
          xi = par$xi, omega = par$omega,
          alpha = par$alpha, solver = "RFB"
        )
      }
    },
    chisq = {
      df <- m2$df
      mir <- isTRUE(m2$mirror)
      function(p) {
        x <- stats::qchisq(p, df = df)
        if (mir) -x else x
      }
    },
    stop("unknown margin2 type")
  )

  # --- sample from Gaussian copula --------------------------------------
  cop <- copula::normalCopula(param = copula_rho, dim = 2)
  u <- draw_u(n, cop, eps) # safe uniforms

  # --- transform ---------------------------------------------------------
  eps_mat <- matrix(NA_real_, n, 2)
  eps_mat[, 1] <- q1(u[, 1])
  eps_mat[, 2] <- q2(u[, 2])
  colnames(eps_mat) <- c("e1", "e2")
  eps_mat
}


# --------------------------------------------------------------------------
# MAIN: simulate_all_conditions_var1
# --------------------------------------------------------------------------
simulate_one_replication <- function(cond_row, rep_i, output_dir, burnin = 30) {
  mu_vec <- c(0, 0) # intercepts 0 (we centred)
  phi_mat <- cond_row$phi_matrix[[1]]
  T_time <- cond_row$T
  T_sim <- T_time + burnin

  # residuals
  eps_mat <- generate_residuals(
    n          = T_sim,
    m1         = cond_row$margin_info[[1]]$margin1,
    m2         = cond_row$margin_info[[1]]$margin2,
    copula_rho = cond_row$rho
  )

  # simulate VAR(1)
  y <- matrix(0, nrow = T_sim, ncol = 2)
  for (t in 2:T_sim) {
    y[t, ] <- mu_vec + phi_mat %*% y[t - 1, ] + eps_mat[t, ]
  }
  y_final <- y[(burnin + 1):T_sim, , drop = FALSE]

  # ---------- package output list ----------------------------------------
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
    rho = cond_row$rho,

    # minimal true‑param bundle for bias calc
    true_params = list(
      mu = mu_vec,
      phi = phi_mat,
      copula_rho = cond_row$rho,
      margin1 = cond_row$margin_info[[1]]$margin1,
      margin2 = cond_row$margin_info[[1]]$margin2
    )
  )

  fn <- file.path(
    output_dir,
    sprintf(
      "sim_data_cond%03d_rep%03d.rds",
      cond_row$condition_id, rep_i
    )
  )
  saveRDS(out, fn)
}

# top‑level factory --------------------------------------------------------
simulate_all_conditions_var1 <- function(sim_conditions_df,
                                         output_dir,
                                         burnin = 30) {
  message("=== Simulating single‑level data ===")
  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]
    msg <- sprintf(
      " Condition %03d  (%s, dir %s, T=%d, rho=%.2f, VAR=%s)",
      cond$condition_id, cond$skew_level,
      cond$direction, cond$T, cond$rho, cond$VARset
    )
    message(msg)
    pb <- utils::txtProgressBar(0, cond$n_reps, style = 3)
    for (r in seq_len(cond$n_reps)) {
      simulate_one_replication(cond, r, output_dir, burnin = burnin)
      utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    cat("\n")
  }
  message("=== Simulation finished ===")
}
