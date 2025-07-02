###########################################################################
# simulate_data.R  –  Multilevel VAR(1) copula data generator  (ASCII only)
# -------------------------------------------------------------------------
# • All parameters are passed in from run_pipeline.R
# • Parallel execution over design rows controlled by the `cores` argument
#   (cores = 1 ⇒ serial).  Cluster workers now receive an explicit export
#   list, so the “could not find function” error cannot occur.
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) stop("Install 'copula'.")
if (!requireNamespace("sn", quietly = TRUE)) stop("Install 'sn'.")
suppressPackageStartupMessages(library(doParallel))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------
# 1.  Helper: skew‑normal parameters with Var = 1
# ------------------------------------------------------------------------
calc_sn_params <- function(alpha, target_var = 1) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(target_var / (1 - 2 * delta^2 / pi))
  xi <- -omega * delta * sqrt(2 / pi)
  list(xi = xi, omega = omega, alpha = alpha)
}

# ------------------------------------------------------------------------
# 2.  Helper: select marginal family & shape from design labels
# ------------------------------------------------------------------------
assign_skew <- function(level, dir_flag) {
  sgn <- function(x) ifelse(x == "+", 1, -1)
  sign1 <- sgn(substr(dir_flag, 1, 1))
  sign2 <- sgn(substr(dir_flag, 2, 2))

  switch(level,
    moderateSN = list(
      margin1 = list(type = "skewnormal", alpha = sign1 * 4),
      margin2 = list(type = "skewnormal", alpha = sign2 * 4)
    ),
    strongSN = list(
      margin1 = list(type = "skewnormal", alpha = sign1 * 9),
      margin2 = list(type = "skewnormal", alpha = sign2 * 9)
    ),
    extremeCHI = list(
      margin1 = list(type = "chisq", df = 1, mirror = (sign1 < 0)),
      margin2 = list(type = "chisq", df = 1, mirror = (sign2 < 0))
    ),
    stop("Unknown skew level: ", level)
  )
}

# ------------------------------------------------------------------------
# 3.  Gaussian‑copula residual generator
# ------------------------------------------------------------------------
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
        df <- m$df
        mir <- isTRUE(m$mirror)
        function(pu) {
          x <- stats::qchisq(pu, df)
          if (mir) -x else x
        }
      },
      stop("Unknown margin type")
    )
  }

  cop <- copula::normalCopula(copula_rho, dim = 2)
  u <- draw_u(n, cop)

  cbind(
    e1 = qfun(m1)(u[, 1]),
    e2 = qfun(m2)(u[, 2])
  )
}

# ------------------------------------------------------------------------
# 4.  Build factorial design grid (reps = # replications per cell)
# ------------------------------------------------------------------------
build_multilevel_design_grid <- function(skew_levels, directions,
                                         T_vals, rhos, var_sets,
                                         reps = 1) {
  expand.grid(
    skew_level = skew_levels,
    direction = directions,
    T = T_vals,
    rho = rhos,
    VARset = names(var_sets),
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    mutate(
      condition_id = row_number(),
      margin_info  = purrr::map2(skew_level, direction, assign_skew),
      phi_matrix   = purrr::map(VARset, \(nm) var_sets[[nm]]),
      n_reps       = reps
    )
}

# ------------------------------------------------------------------------
# 5.  Subject‑level simulator (returns matrix T × 2)
# ------------------------------------------------------------------------
simulate_subject_series <- function(T_time, mu_vec, phi_mat,
                                    margin_info, copula_rho, burnin) {
  eps <- generate_residuals(
    n          = T_time + burnin,
    m1         = margin_info$margin1,
    m2         = margin_info$margin2,
    copula_rho = copula_rho
  )

  y <- matrix(0, T_time + burnin, 2)
  for (tt in 2:nrow(y)) {
    y[tt, ] <- mu_vec + phi_mat %*% y[tt - 1, ] + eps[tt, ]
  }

  y[(burnin + 1):(burnin + T_time), , drop = FALSE]
}

# ------------------------------------------------------------------------
# 6.  Generate ONE replication file for a design cell
# ------------------------------------------------------------------------
simulate_one_replication_ml <- function(cond_row,
                                        N, rep_i, output_dir,
                                        sd_AR, sd_CL, burnin) {
  T_time <- cond_row$T
  phi_fixed <- cond_row$phi_matrix[[1]]

  dat_list <- vector("list", N)
  phi_draws <- array(NA_real_, dim = c(N, 2, 2))

  for (i in seq_len(N)) {
    phi_i <- phi_fixed
    phi_i[1, 1] <- phi_fixed[1, 1] + rnorm(1, 0, sd_AR)
    phi_i[2, 2] <- phi_fixed[2, 2] + rnorm(1, 0, sd_AR)
    phi_i[1, 2] <- phi_fixed[1, 2] + rnorm(1, 0, sd_CL)
    phi_i[2, 1] <- phi_fixed[2, 1] + rnorm(1, 0, sd_CL)

    y_i <- simulate_subject_series(
      T_time, c(0, 0), phi_i,
      cond_row$margin_info[[1]],
      cond_row$rho, burnin
    )

    dat_list[[i]] <- tibble::tibble(
      i  = i,
      t  = seq_len(T_time),
      y1 = y_i[, 1],
      y2 = y_i[, 2]
    )
    phi_draws[i, , ] <- phi_i
  }

  saveRDS(
    list(
      condition_id = cond_row$condition_id,
      rep_i        = rep_i,
      N            = N,
      T            = T_time,
      data         = dplyr::bind_rows(dat_list),
      phi_fixed    = phi_fixed,
      phi_subject  = phi_draws,
      rho          = cond_row$rho,
      margin_info  = cond_row$margin_info[[1]],
      sd_AR        = sd_AR,
      sd_CL        = sd_CL,
      burnin       = burnin
    ),
    file.path(
      output_dir,
      sprintf(
        "sim_data_ml_cond%03d_rep%03d.rds",
        cond_row$condition_id, rep_i
      )
    )
  )
}

# ------------------------------------------------------------------------
# 7.  Wrapper: all conditions, optional parallel execution
# ------------------------------------------------------------------------
simulate_all_conditions_ml <- function(sim_conditions_df,
                                       output_dir,
                                       N_per_subject,
                                       burnin,
                                       sd_AR, sd_CL,
                                       cores = 1) {
  dir.create(output_dir, FALSE, TRUE)

  run_condition <- function(cond) {
    for (rep_i in seq_len(cond$n_reps)) {
      message(sprintf(
        "Cond %03d Rep %02d | %s%s  T=%d  rho=%.2f  VAR=%s  N=%d",
        cond$condition_id, rep_i,
        cond$skew_level, cond$direction,
        cond$T, cond$rho, cond$VARset, N_per_subject
      ))
      simulate_one_replication_ml(
        cond_row   = cond,
        N          = N_per_subject,
        rep_i      = rep_i,
        output_dir = output_dir,
        sd_AR      = sd_AR,
        sd_CL      = sd_CL,
        burnin     = burnin
      )
    }
    NULL
  }

  if (cores <= 1) {
    lapply(
      seq_len(nrow(sim_conditions_df)),
      \(idx) run_condition(sim_conditions_df[idx, ])
    )
  } else {
    cl <- parallel::makeCluster(cores)
    registerDoParallel(cl)
    on.exit({
      stopCluster(cl)
      registerDoSEQ()
    })

    foreach::foreach(
      idx = seq_len(nrow(sim_conditions_df)),
      .packages = c("dplyr", "tibble", "copula", "sn"),
      ## explicit exports so workers can see them
      .export = c(
        "run_condition",
        "simulate_one_replication_ml",
        "simulate_subject_series",
        "generate_residuals",
        "calc_sn_params",
        "assign_skew",
        "N_per_subject", "burnin",
        "sd_AR", "sd_CL", "output_dir"
      )
    ) %dopar% {
      run_condition(sim_conditions_df[idx, ])
    }
  }

  message(">>> Data generation complete.")
}
