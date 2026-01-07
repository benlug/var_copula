###########################################################################
# simulate_data_SEM.R — Study 5 (SEM): Simulation
#
# Purpose
#   Generate datasets for Study 5, which contrasts two placements of
#   signed/shifted Exponential margins:
#     A_indicator: Exponential margins at the *measurement* layer
#     B_latent   : Exponential margins at the *innovation* layer
#
# Design goals
#   - Deterministic, per-(condition_id, rep_id) seeding (resume-safe)
#   - Skip existing simulation files unless overwrite = TRUE
#   - Correct copula behaviour under sign-flips ("+-" direction)
#
# Critical technical note on sign flips
#   If we want a left-skewed signed Exponential residual via reflection,
#   we must *not* use e = 1 - Qexp(u) with the same copula-uniform u.
#   That transformation is monotone decreasing in u and rotates the copula
#   (for Gaussian copula it flips the sign of dependence when only one
#   margin is flipped). To preserve the copula parameter rho, the correct
#   increasing quantile is:
#       e(u) = 1 - Qexp(1-u)
#
#   This mirrors the fix used for the mirrored chi-square condition in
#   simulate_data.R (Studies 1–3).
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

# ------------------------------------------------------------------------
# 1) Gaussian copula uniforms U in (0,1)
# ------------------------------------------------------------------------
r_gauss_copula_u <- function(n, rho, eps = 1e-12) {
  # MASS is a recommended package; use namespace call to avoid attaching.
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for MASS::mvrnorm().")
  }

  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  Z <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  U <- cbind(stats::pnorm(Z[, 1]), stats::pnorm(Z[, 2]))

  # Clamp away from {0,1} to avoid extreme quantiles.
  U <- pmin(pmax(U, eps), 1 - eps)
  matrix(as.numeric(U), nrow = n, ncol = 2, byrow = FALSE)
}

# ------------------------------------------------------------------------
# 2) Signed, standardized Exponential quantiles (rate = 1)
# ------------------------------------------------------------------------
# Right-skew (dir = +1):  e = X - 1,  X ~ Exp(1), support e >= -1
# Left-skew  (dir = -1):  e = 1 - X,  X ~ Exp(1), support e <=  1
#
# IMPORTANT: to preserve the copula parameter rho under left-skew, we must
# generate left-skew via the *increasing* quantile Q_{1-X}(u) = 1-Q_X(1-u).
q_exp_signed <- function(u_vec, dir, eps = 1e-12) {
  u_vec <- pmin(pmax(u_vec, eps), 1 - eps)
  if (dir > 0L) {
    stats::qexp(u_vec, rate = 1) - 1
  } else {
    1 - stats::qexp(1 - u_vec, rate = 1)
  }
}

# ------------------------------------------------------------------------
# 3) One dataset
# ------------------------------------------------------------------------
simulate_one <- function(sem_study, T, B, rho, skew_dir) {
  mu <- c(0, 0)
  skew_dir <- as.integer(skew_dir) # length 2: c(+1/-1, +1/-1)

  if (identical(sem_study, "A_indicator")) {
    # Study A: Gaussian latent state; Exponential-skew measurement error
    state <- matrix(NA_real_, nrow = T, ncol = 2)
    z_state <- matrix(NA_real_, nrow = T, ncol = 2)

    # Align with the Stan model: x_1 = mu + z_1.
    s <- mu + stats::rnorm(2, 0, 1)
    state[1, ] <- s
    z_state[1, ] <- s - mu

    if (T > 1) {
      for (t in 2:T) {
        zt <- stats::rnorm(2, 0, 1)
        s <- mu + as.vector(B %*% s) + zt
        state[t, ] <- s
        z_state[t, ] <- zt
      }
    }
    stopifnot(all(dim(state) == c(T, 2)))

    # Copula uniforms for measurement layer
    U <- r_gauss_copula_u(T, rho)
    e <- cbind(
      q_exp_signed(U[, 1], skew_dir[1]),
      q_exp_signed(U[, 2], skew_dir[2])
    )
    e <- matrix(as.numeric(e), nrow = T, ncol = 2, byrow = FALSE)
    stopifnot(all(dim(e) == c(T, 2)))

    y <- state + e
    stopifnot(all(dim(y) == c(T, 2)))

    return(list(
      data       = data.frame(y1 = y[, 1], y2 = y[, 2]),
      T          = T,
      rho        = rho,
      sem_study  = "A_indicator",
      phi_matrix = B,
      skew_signs = skew_dir,
      truth      = list(state = state, meas_error = e, innovations = z_state)
    ))
  }

  if (identical(sem_study, "B_latent")) {
    # Study B: VAR(1) with Exponential-skew innovations (copula at innovation layer)
    y <- matrix(0, nrow = T, ncol = 2)
    y[1, ] <- c(0, 0) # fixed initial condition (matches Stan model)

    U <- r_gauss_copula_u(T - 1, rho)
    z <- cbind(
      q_exp_signed(U[, 1], skew_dir[1]),
      q_exp_signed(U[, 2], skew_dir[2])
    )
    z <- matrix(as.numeric(z), nrow = T - 1, ncol = 2, byrow = FALSE)
    stopifnot(all(dim(z) == c(T - 1, 2)))

    for (t in 2:T) {
      y[t, ] <- mu + as.vector(B %*% y[t - 1, ]) + z[t - 1, ]
    }
    y <- matrix(as.numeric(y), nrow = T, ncol = 2, byrow = FALSE)
    stopifnot(all(dim(y) == c(T, 2)))

    innovations <- rbind(c(NA_real_, NA_real_), z)

    return(list(
      data       = data.frame(y1 = y[, 1], y2 = y[, 2]),
      T          = T,
      rho        = rho,
      sem_study  = "B_latent",
      phi_matrix = B,
      skew_signs = skew_dir,
      truth      = list(state = NULL, meas_error = NULL, innovations = innovations)
    ))
  }

  stop("Unknown sem_study: ", sem_study)
}

# ------------------------------------------------------------------------
# 4) Driver over all conditions (resume-aware)
# ------------------------------------------------------------------------
simulate_all_conditions_SEM <- function(sim_conditions_df,
                                        output_dir,
                                        start_condition = 1,
                                        start_rep = 1,
                                        overwrite = FALSE,
                                        seed_base = 2025000L,
                                        verbose = TRUE) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose) {
    message("=== Simulating Study 5 SEM datasets ===")
    message("Output dir: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
  }

  # Safety: ensure deterministic RNG kind across R versions (optional).
  # This mirrors common reproducibility practice and is harmless if left as-is.
  suppressWarnings({
    RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
  })

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]
    cid <- as.integer(cond$condition_id)

    if (cid < start_condition) next
    first_rep <- if (cid == start_condition) start_rep else 1

    if (verbose) {
      msg <- sprintf(
        " Condition %03d  (sem=%s dir=%s T=%d rho=%.2f)",
        cid, cond$sem_study, cond$direction, cond$T, cond$rho
      )
      message(msg)
    }

    pb <- utils::txtProgressBar(min = 1, max = cond$n_reps,
                                initial = first_rep, style = 3)

    for (r in seq(from = first_rep, to = cond$n_reps)) {
      target <- file.path(output_dir, sprintf("sim_data_cond%03d_rep%03d.rds", cid, r))
      if (!overwrite && file.exists(target)) {
        utils::setTxtProgressBar(pb, r)
        next
      }

      # Deterministic per-(condition, replication) seeding
      seed <- as.integer(seed_base + cid * 10000L + as.integer(r))
      if (is.na(seed) || seed <= 0L) seed <- as.integer(abs(seed) + 1L)
      set.seed(seed)

      sim <- simulate_one(
        sem_study = cond$sem_study,
        T         = cond$T,
        B         = cond$phi_matrix[[1]],
        rho       = cond$rho,
        skew_dir  = cond$skew_dir[[1]]
      )

      # Add provenance
      sim$condition_id <- cid
      sim$rep_id <- as.integer(r)
      sim$seed <- seed
      sim$direction <- cond$direction

      saveRDS(sim, target)
      utils::setTxtProgressBar(pb, r)
    }

    close(pb)
    cat("\n")
  }

  if (verbose) message("=== Simulation finished ===")
}
