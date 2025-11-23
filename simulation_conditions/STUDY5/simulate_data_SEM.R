# simulate_data_SEM.R — robust shapes for SEM Study 5

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

# Draw U ~ Gaussian copula (rho) → strict n×2 matrix in (0,1)
r_gauss_copula_u <- function(n, rho) {
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  Z <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma) # do not attach MASS
  U <- cbind(pnorm(Z[, 1]), pnorm(Z[, 2]))
  # force strict matrix shape, even if n == 1
  matrix(as.numeric(U), nrow = n, ncol = 2, byrow = FALSE)
}

# Signed standardized exponential (rate=1)
# dir = +1 (right): e = x - 1  (support e ≥ -1)
# dir = -1 (left) : e = 1 - x  (support e ≤  1)
q_exp_signed <- function(u_vec, dir) {
  x <- stats::qexp(u_vec, rate = 1)
  ifelse(dir > 0, x - 1, 1 - x)
}

simulate_one <- function(sem_study, T, B, rho, skew_dir) {
  mu <- c(0, 0)
  skew_dir <- as.integer(skew_dir) # length 2: c(+1/-1, +1/-1)

  if (identical(sem_study, "A_indicator")) {
    # ---------------- Study A: Gaussian state; exponential-skew measurement error
    state <- matrix(NA_real_, nrow = T, ncol = 2)
    s <- mu
    for (t in seq_len(T)) {
      if (t > 1) {
        # B %*% s returns 2×1; coerce to vector before + to avoid matrix carry-over
        s <- mu + as.vector(B %*% s) + rnorm(2, 0, 1)
      }
      state[t, ] <- s
    }
    stopifnot(all(dim(state) == c(T, 2)))

    # Copula uniforms for measurement layer (NO mirroring; use sign in quantile)
    U <- r_gauss_copula_u(T, rho)
    e <- cbind(
      q_exp_signed(U[, 1], skew_dir[1]),
      q_exp_signed(U[, 2], skew_dir[2])
    )
    e <- matrix(as.numeric(e), nrow = T, ncol = 2, byrow = FALSE)
    stopifnot(all(dim(e) == c(T, 2)))

    # Now both are guaranteed T×2 numeric
    y <- matrix(as.numeric(state), nrow = T, ncol = 2) +
      matrix(as.numeric(e), nrow = T, ncol = 2)
    stopifnot(all(dim(y) == c(T, 2)))

    return(list(
      data        = data.frame(y1 = y[, 1], y2 = y[, 2]),
      T           = T,
      rho         = rho,
      sem_study   = "A_indicator",
      phi_matrix  = B,
      skew_signs  = skew_dir,
      truth       = list(state = state, meas_error = e, innovations = NULL)
    ))
  }

  if (identical(sem_study, "B_latent")) {
    # ---------------- Study B: VAR(1) with exponential-skew innovations (copula at innovation layer)
    y <- matrix(0, nrow = T, ncol = 2)
    y[1, ] <- c(0, 0) # η1 = 0 per study spec

    # Copula uniforms for innovation layer (NO mirroring; use sign in quantile)
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
      data        = data.frame(y1 = y[, 1], y2 = y[, 2]),
      T           = T,
      rho         = rho,
      sem_study   = "B_latent",
      phi_matrix  = B,
      skew_signs  = skew_dir,
      truth       = list(state = NULL, meas_error = NULL, innovations = innovations)
    ))
  }

  stop("Unknown sem_study: ", sem_study)
}

simulate_all_conditions_SEM <- function(sim_conditions_df, output_dir,
                                        start_condition = 1, start_rep = 1) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (i in seq_len(nrow(sim_conditions_df))) {
    if (sim_conditions_df$condition_id[i] < start_condition) next

    sem_study <- sim_conditions_df$sem_study[i]
    T <- sim_conditions_df$T[i]
    rho <- sim_conditions_df$rho[i]
    B <- sim_conditions_df$phi_matrix[[i]]
    skew_dir <- sim_conditions_df$skew_dir[[i]]
    n_reps <- sim_conditions_df$n_reps[i]

    for (r in seq_len(n_reps)) {
      if (sim_conditions_df$condition_id[i] == start_condition && r < start_rep) next

      set.seed(1e6 + 1000L * i + r)
      sim <- simulate_one(sem_study, T, B, rho, skew_dir)

      fn <- file.path(output_dir, sprintf(
        "sim_data_cond%03d_rep%03d.rds",
        sim_conditions_df$condition_id[i], r
      ))
      saveRDS(sim, fn)
    }
  }
  message("=== Simulation finished ===")
}
