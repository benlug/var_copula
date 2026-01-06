# simulate_data_multilevel.R — EG DGP, minimal pooling (random intercepts only)
if (!requireNamespace("copula", quietly = TRUE)) stop("Install package 'copula'.")

# reuse your generator for Exponential margins + Gaussian copula
generate_residuals <- function(n, m1, m2, copula_rho, eps = 1e-9) {
  draw_u <- function(n, cop) {
    u <- copula::rCopula(n, cop)
    while (any(u <= eps | u >= 1 - eps)) {
      bad <- unique(which(u <= eps | u >= 1 - eps, arr.ind = TRUE)[, 1])
      u[bad, ] <- copula::rCopula(length(bad), cop)
    }
    u
  }
  qfun <- function(m) {
    switch(m$type,
      exponential = {
        rate <- m$rate
        mean_exp <- 1 / rate
        sd_exp <- 1 / rate
        mir <- isTRUE(m$mirror)
        function(pu) {
          x_raw <- stats::qexp(pu, rate)
          x_std <- (x_raw - mean_exp) / sd_exp
          if (mir) -x_std else x_std
        }
      },
      stop("Only 'exponential' is supported in the ML EG DGP.")
    )
  }
  cop <- copula::normalCopula(copula_rho, dim = 2)
  u <- draw_u(n, cop)
  cbind(e1 = qfun(m1)(u[, 1]), e2 = qfun(m2)(u[, 2]))
}

# simple projection to ensure stationarity if needed
.project_if_needed <- function(Phi, alpha = 0.995) {
  ev <- eigen(Phi, only.values = TRUE)$values
  r <- max(Mod(ev))
  if (is.na(r) || r < 1) {
    return(Phi)
  }
  Phi * (alpha / r)
}

simulate_one_replication_ml <- function(cond_row, rep_i, output_dir) {
  # Deterministic per-(condition, replication) seeding.
  #
  # Why this matters:
  #   - If you resume a run with START_COND/START_REP while skipping
  #     already-written .rds files, using only a single global set.seed()
  #     can silently generate repeated datasets under different rep IDs.
  #   - Seeding inside each replication makes simulation reproducible and
  #     resume-safe.
  seed_sim <- 700000L + as.integer(cond_row$condition_id) * 1000L + as.integer(rep_i)
  set.seed(seed_sim)
  N <- cond_row$N_units
  Tn <- cond_row$T
  burnin <- cond_row$burnin
  T_sim <- Tn + burnin

  Phi <- .project_if_needed(cond_row$phi_matrix[[1]])
  rho <- cond_row$rho
  tau <- cond_row$tau_mu[[1]] # length 2
  stopifnot(length(tau) == 2)

  # random intercepts (minimal pooling DGP)
  mu_mat <- cbind(
    rnorm(N, 0, tau[1]),
    rnorm(N, 0, tau[2])
  )

  # innovations per unit (same margin/cor across units)
  m1 <- cond_row$margin_info[[1]]$margin1
  m2 <- cond_row$margin_info[[1]]$margin2

  all_rows <- vector("list", N)
  for (i in seq_len(N)) {
    eps <- generate_residuals(T_sim, m1, m2, rho)
    y <- matrix(0, T_sim, 2)
    y[1, ] <- eps[1, ]
    for (t in 2:T_sim) {
      y[t, ] <- mu_mat[i, ] + Phi %*% y[t - 1, ] + eps[t, ]
    }
    y_fin <- y[(burnin + 1):T_sim, , drop = FALSE]
    all_rows[[i]] <- data.frame(
      i = i, t = seq_len(Tn),
      y1 = y_fin[, 1], y2 = y_fin[, 2]
    )
  }
  panel <- dplyr::bind_rows(all_rows)

  # helper to coerce anything to a single TRUE/FALSE
  .flag <- function(x) isTRUE(as.logical(x)[1])

  out <- list(
    condition_id = cond_row$condition_id,
    rep_i = rep_i,
    seed_sim = seed_sim,
    N = N, T = Tn,
    data = panel,
    phi_matrix = Phi,
    rho = rho,
    true_params = list(
      mu_bar = c(0, 0),
      tau_mu = tau,
      mu_mat = mu_mat,
      phi = Phi,
      copula_rho = rho,
      margin1 = m1,
      margin2 = m2,
      sigma_exp = c(1, 1),
      skew_dir = c(
        if (.flag(m1$mirror)) -1 else 1,
        if (.flag(m2$mirror)) -1 else 1
      )
    )
  )

  saveRDS(out, file.path(
    output_dir, sprintf(
      "sim_dataML_cond%03d_rep%03d.rds",
      cond_row$condition_id, rep_i
    )
  ))
}

simulate_all_conditions_var1_ml <- function(sim_conditions_df,
                                            output_dir,
                                            start_condition = 1,
                                            start_rep = 1,
                                            overwrite = FALSE) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating MULTILEVEL (minimal pooling) data ===")

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]
    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) start_rep else 1

    message(sprintf(
      " Condition %03d  (%s dir=%s N=%d T=%d ρ=%.2f VAR=%s τμ=(%.2f,%.2f))",
      cond$condition_id, cond$dgp_level, cond$direction,
      cond$N_units, cond$T, cond$rho, cond$VARset,
      cond$tau_mu[[1]][1], cond$tau_mu[[1]][2]
    ))

    pb <- utils::txtProgressBar(
      min = 1, max = cond$n_reps,
      initial = first_rep, style = 3
    )
    for (r in seq(from = first_rep, to = cond$n_reps)) {
      f <- file.path(
        output_dir,
        sprintf(
          "sim_dataML_cond%03d_rep%03d.rds",
          cond$condition_id, r
        )
      )
      if (!overwrite && file.exists(f)) {
        utils::setTxtProgressBar(pb, r)
        next
      }
      simulate_one_replication_ml(cond, r, output_dir)
      utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    cat("\n")
  }
  message("=== ML simulation finished ===")
}
