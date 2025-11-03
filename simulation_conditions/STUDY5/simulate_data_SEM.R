###########################################################################
# simulate_data_SEM.R — SEM Study A/B with Exponential margins
# Study A: indicator-skew (ε ~ Exp^± with copula ρ, ζ ~ N)
# Study B: latent-skew    (ζ ~ Exp^± with copula ρ, ε ≡ 0)
# Exponential quantiles follow Study 4 note. :contentReference[oaicite:3]{index=3}
###########################################################################

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Install package 'copula'.")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Exponential(+): E(+) = X - 1, X~Exp(1). Left-skew via mirror E(-) = -E(+).
# Quantiles (standardized, right/left):
q_exp <- function(u, sign = +1L) {
  stopifnot(all(u > 0 & u < 1))
  if (sign == +1L) { # right-skew
    -1 - log1p(-u)
  } else { # left-skew
    1 + log(u)
  }
}

# draw 2D active-layer residuals via Gaussian copula then marginal Q(+/−)
draw_active <- function(n, rho, signs) {
  cop <- copula::normalCopula(rho, dim = 2)
  U <- copula::rCopula(n, cop)
  cbind(
    q_exp(U[, 1], signs[1]),
    q_exp(U[, 2], signs[2])
  )
}

simulate_one_SEM <- function(cond_row, rep_i, output_dir) {
  Tt <- cond_row$T
  B <- cond_row$phi_matrix[[1]]
  rho_act <- cond_row$rho
  sgn <- unlist(cond_row$skew_dir)

  # State innovations and measurement errors per study
  if (cond_row$sem_study == "A_indicator") {
    # Study A: ζ ~ N(0, I); ε ~ Exp^± (copula ρ) at measurement layer.
    z <- matrix(rnorm(Tt * 2, 0, 1), Tt, 2)
    eta <- matrix(0, Tt, 2)
    for (t in 2:Tt) eta[t, ] <- as.numeric(B %*% eta[t - 1, ]) + z[t, ]
    eps <- draw_active(Tt, rho_act, sgn)
    y <- eta + eps
    active_layer <- "measurement"
  } else {
    # Study B: ζ ~ Exp^± (copula ρ) at state layer; ε ≡ 0; y = η deterministically.
    z <- draw_active(Tt, rho_act, sgn)
    eta <- matrix(0, Tt, 2)
    for (t in 2:Tt) eta[t, ] <- as.numeric(B %*% eta[t - 1, ]) + z[t, ]
    y <- eta
    active_layer <- "state"
  }

  out <- list(
    condition_id = cond_row$condition_id,
    rep_i = rep_i,
    T = Tt,
    data = data.frame(t = seq_len(Tt), y1 = y[, 1], y2 = y[, 2]),
    phi_matrix = B,
    rho = rho_act,
    sem_study = cond_row$sem_study,
    direction = cond_row$direction,
    active_layer = active_layer,
    skew_signs = sgn,
    true_params = list(
      mu = c(0, 0),
      phi = B,
      rho = rho_act,
      sigma_exp = c(1, 1), # standardized margins
      sem_study = cond_row$sem_study
    )
  )

  saveRDS(
    out,
    file.path(output_dir, sprintf(
      "sim_data_cond%03d_rep%03d.rds",
      cond_row$condition_id, rep_i
    ))
  )
}

simulate_all_conditions_SEM <- function(sim_conditions_df, output_dir,
                                        start_condition = 1, start_rep = 1,
                                        overwrite = FALSE) {
  dir.create(output_dir, FALSE, TRUE)
  message("=== Simulating SEM data (short design) ===")

  for (rr in seq_len(nrow(sim_conditions_df))) {
    cond <- sim_conditions_df[rr, ]

    if (cond$condition_id < start_condition) next
    first_rep <- if (cond$condition_id == start_condition) start_rep else 1L

    message(sprintf(
      " Condition %03d  (%s dir=%s T=%d rho=%.2f)",
      cond$condition_id, cond$sem_study, cond$direction,
      cond$T, cond$rho
    ))
    pb <- utils::txtProgressBar(min = 1, max = cond$n_reps, initial = first_rep, style = 3)

    for (r in seq(from = first_rep, to = cond$n_reps)) {
      fp <- file.path(output_dir, sprintf(
        "sim_data_cond%03d_rep%03d.rds",
        cond$condition_id, r
      ))
      if (!overwrite && file.exists(fp)) {
        utils::setTxtProgressBar(pb, r)
        next
      }
      simulate_one_SEM(cond, r, output_dir)
      utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    cat("\n")
  }
  message("=== Simulation finished ===")
}
