###########################################################################
# simulate_data.R – Study III: Normal margins + Clayton copula DGP
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Install package 'copula'.")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----- quantile functions for margins (standard normal) ------------------
qfun_for_margin <- function(m) {
  switch(m$type,
    normal = function(pu) stats::qnorm(pu),
    stop("unknown margin type: ", m$type)
  )
}

# ----- robust sampler for copula uniforms --------------------------------
draw_u <- function(n, cop) {
  eps <- 1e-9
  u <- copula::rCopula(n, cop)
  while (any(u <= eps | u >= 1 - eps)) {
    bad <- unique(which(u <= eps | u >= 1 - eps, arr.ind = TRUE)[, 1])
    u[bad, ] <- copula::rCopula(length(bad), cop)
  }
  u
}

generate_residuals <- function(n, m1, m2, clayton_theta) {
  cop <- copula::claytonCopula(clayton_theta, dim = 2)
  u <- draw_u(n, cop)
  q1 <- qfun_for_margin(m1)
  q2 <- qfun_for_margin(m2)
  cbind(e1 = q1(u[, 1]), e2 = q2(u[, 2]))
}

simulate_one_replication <- function(cond_row, rep_i, output_dir, burnin = 30) {
  mu_vec <- c(0, 0)
  phi_mat <- cond_row$phi_matrix[[1]]
  T_time <- cond_row$T
  T_sim <- T_time + burnin

  eps_mat <- generate_residuals(
    n             = T_sim,
    m1            = cond_row$margin_info[[1]]$margin1,
    m2            = cond_row$margin_info[[1]]$margin2,
    clayton_theta = cond_row$theta
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
    data = data.frame(i = 1, t = seq_len(T_time), y1 = y_final[, 1], y2 = y_final[, 2]),
    phi_matrix = phi_mat,
    theta = cond_row$theta, # convenience copy
    true_params = list(
      mu           = mu_vec,
      phi          = phi_mat,
      copula_type  = "clayton",
      copula_theta = cond_row$theta,
      margin1      = cond_row$margin_info[[1]]$margin1,
      margin2      = cond_row$margin_info[[1]]$margin2
    )
  )

  saveRDS(
    out,
    file.path(output_dir, sprintf("sim_data_cond%03d_rep%03d.rds", cond_row$condition_id, rep_i))
  )
}

simulate_all_conditions_var1 <- function(sim_conditions_df, output_dir, burnin = 30,
                                         start_condition = 1, start_rep = 1,
                                         overwrite = FALSE) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating Study III data ===")

  # Harmonize for logging
  if ("dgp_level" %in% names(sim_conditions_df) && !"skew_level" %in% names(sim_conditions_df)) {
    names(sim_conditions_df)[names(sim_conditions_df) == "dgp_level"] <- "skew_level"
  }

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]
    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) start_rep else 1

    message(sprintf(
      " Condition %03d (%s T=%d θ=%.2f VAR=%s)",
      cond$condition_id, cond$skew_level, cond$T, cond$theta, cond$VARset
    ))
    pb <- utils::txtProgressBar(min = 1, max = cond$n_reps, initial = first_rep, style = 3)

    for (r in seq(from = first_rep, to = cond$n_reps)) {
      file_target <- file.path(output_dir, sprintf("sim_data_cond%03d_rep%03d.rds", cond$condition_id, r))
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
