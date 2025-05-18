###########################################################################
# simulate_data.R (Generalized Version)
#
# Simulates multilevel VAR(1) data based on conditions provided.
# Innovations can be left/right Skew-Normal, zero-mean.
# Includes options for random effects (handled via SDs).
# Called by run_pipeline.R
###########################################################################

# --- Libraries ---
# Ensure necessary libraries are loaded by the calling script (run_pipeline.R)
# Requires: MASS, dplyr, purrr (for pipeline), sn

# --- Helper Functions ---

#' Draw from Gaussian Copula
#' @param n Number of samples
#' @param copula_info List containing copula type ("gaussian") and params (list(rho=...))
#' @return Matrix of uniforms (n x 2)
draw_copula <- function(n = 1, copula_info) {
    # Ensure MASS is loaded for mvrnorm
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' needed for mvrnorm.")
    # Ensure base stats package is available for pnorm
    if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' needed for pnorm.")

    if (copula_info$copula == "gaussian") {
        rho <- copula_info$params$rho
        if (is.null(rho) || !is.numeric(rho) || length(rho) != 1 || abs(rho) >= 1) {
            stop("Correlation rho must be a single numeric value between -1 and 1.")
        }
        Sigma_z <- matrix(c(1, rho, rho, 1), nrow = 2)
        mu_z <- c(0, 0)

        Z <- tryCatch(MASS::mvrnorm(n, mu = mu_z, Sigma = Sigma_z),
                      error = function(e) stop("mvrnorm failed in draw_copula: ", e$message))
        # Ensure output is always a matrix
        if (n == 1 && is.vector(Z)) Z <- matrix(Z, nrow = 1)

        U <- stats::pnorm(Z)
        return(U)
    } else {
        stop("Unsupported copula type: ", copula_info$copula)
    }
}


#' Draw from Specified Marginal Distribution using Inverse CDF
#' @param u Uniform value(s) (0, 1)
#' @param dist_params List defining the distribution (dist="skewnormal", params=list(xi=.., omega=.., alpha=..))
#' @return Value(s) drawn from the specified marginal distribution
draw_from_margin <- function(u, dist_params) {
    # Ensure sn package is loaded for qsn
    if (!requireNamespace("sn", quietly = TRUE)) stop("Package 'sn' needed for qsn.")

    if (dist_params$dist != "skewnormal") {
        stop("This simulation currently only supports 'skewnormal' margins.")
    }
    sim_params <- dist_params$params
    xi_val <- sim_params$xi
    omega_val <- sim_params$omega
    alpha_val <- sim_params$alpha

    # Basic validation
    if (is.null(xi_val) || is.null(omega_val) || is.null(alpha_val)) stop("Missing skew-normal parameters (xi, omega, alpha).")
    if (!is.numeric(xi_val) || !is.numeric(omega_val) || !is.numeric(alpha_val)) stop("Skew-normal parameters must be numeric.")
    if (omega_val <= 0) stop("Skew-Normal scale (omega) must be positive.")

    # Clamp uniform values to avoid issues at boundaries with qsn
    eps <- 1e-9 # Small epsilon to avoid exact 0 or 1
    u_clamp <- pmax(eps, pmin(1 - eps, u)) # Vectorized clamp

    # Use qsn from 'sn' package
    vals <- tryCatch(
                sn::qsn(u_clamp, xi = xi_val, omega = omega_val, alpha = alpha_val, solver = "RFB"),
                error = function(e) {
                    # Add more context to the error if it still fails
                    stop(paste0("qsn failed in draw_from_margin (using solver='RFB') with inputs: ",
                                "u_clamp=", round(u_clamp, 5), ", xi=", round(xi_val, 3),
                                ", omega=", round(omega_val, 3), ", alpha=", round(alpha_val, 3),
                                ". Original error: ", e$message))
                }
            )
    return(vals)
}


# --- Main Simulation Function ---

#' Simulate data for all conditions specified in a data frame
#'
#' @param sim_conditions_df Data frame where each row defines a simulation condition.
#'                          Must contain columns like N, T, margin1, margin2, copula_choice,
#'                          phi_base elements, sigma_mu_vec, sigma_phi_vec, n_reps, etc.
#'                          and a unique 'condition_id'. `has_random_effects` should also be present.
#' @param output_dir Directory where the simulation RDS files will be saved.
#' @return Invisibly returns NULL. Saves files to disk.
simulate_all_conditions <- function(sim_conditions_df, output_dir) {

    # Ensure necessary packages for internal operations are available
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' needed.")

    if (!dir.exists(output_dir)) {
        warning("Output directory did not exist, creating: ", output_dir)
        tryCatch(dir.create(output_dir, recursive = TRUE),
                  error = function(e) stop("Failed to create output directory: ", e$message))
    }

    total_conditions <- nrow(sim_conditions_df)
    cat(sprintf("Starting simulation for %d conditions...\n", total_conditions))

    # Loop through each condition (row) in the data frame
    for (cond_idx in seq_len(total_conditions)) {
        cond <- sim_conditions_df[cond_idx, ]
        current_cond_id <- cond$condition_id # Use the unique ID

        cat(sprintf("\n--- Starting Condition ID %d/%d ---\n", current_cond_id, total_conditions))
        # Print key parameters for this condition
        # Safely access nested list elements
        alpha_val <- tryCatch(cond$margin1[[1]]$params$alpha, error = function(e) NA)
        cat(sprintf("  N=%d, T=%d, Alpha=%.1f, AR=%.1f, Cross=%.1f, RE=%s\n",
                    cond$N, cond$T, alpha_val, cond$phi11_base,
                    cond$phi12_base, cond$has_random_effects))

        # Extract Scenario Parameters for this row
        N_val <- cond$N; T_val <- cond$T
        mu_global_re_mean <- c(cond$mu_global_1, cond$mu_global_2)
        phi11_base <- cond$phi11_base; phi12_base <- cond$phi12_base
        phi21_base <- cond$phi21_base; phi22_base <- cond$phi22_base
        # sigma_mu_vec and sigma_phi_vec are lists, access with [[1]]
        sigma_mu_vec <- cond$sigma_mu_vec[[1]]
        sigma_phi_vec <- cond$sigma_phi_vec[[1]]
        n_reps <- cond$n_reps
        margin1_info <- cond$margin1[[1]]
        margin2_info <- cond$margin2[[1]]
        copula_info <- cond$copula_choice[[1]]
        burnin <- cond$burnin

        # Validate parameters
        if (length(sigma_mu_vec) != 2) stop(sprintf("Condition %d: sigma_mu_vec must have length 2.", current_cond_id))
        if (length(sigma_phi_vec) != 4) stop(sprintf("Condition %d: sigma_phi_vec must have length 4.", current_cond_id))
        if (burnin < 0) stop(sprintf("Condition %d: Burn-in must be non-negative.", current_cond_id))
        if (!is.list(margin1_info) || !is.list(margin2_info) || !is.list(copula_info)) {
             stop(sprintf("Condition %d: Margin/Copula info is not correctly formatted as lists.", current_cond_id))
        }

        T_total <- T_val + burnin

        # --- Loop through replications for the current condition ---
        for (rep_i in seq_len(n_reps)) {
            # Generate a unique seed for each rep within each condition
            data_seed <- 10000 + current_cond_id * 1000 + rep_i # Ensure uniqueness
            set.seed(data_seed)
            # Less verbose output within replication loop
            if (rep_i == 1 || rep_i %% 50 == 0 || rep_i == n_reps) {
                cat(sprintf("  Replication %d/%d (Cond ID: %d, Seed: %d)...\n",
                            rep_i, n_reps, current_cond_id, data_seed))
            }

            # --- Generate Random Effects (will be zero if SDs are zero) ---
            mu_i_raw <- matrix(stats::rnorm(N_val * 2), ncol = 2)
            dev_phi_raw <- matrix(stats::rnorm(N_val * 4), ncol = 4)

            # Scale raw deviates by the SDs
            mu_i <- matrix(0, nrow = N_val, ncol = 2)
            mu_i[, 1] <- mu_global_re_mean[1] + mu_i_raw[, 1] * sigma_mu_vec[1]
            mu_i[, 2] <- mu_global_re_mean[2] + mu_i_raw[, 2] * sigma_mu_vec[2]

            dev_phi11 <- dev_phi_raw[, 1] * sigma_phi_vec[1]
            dev_phi12 <- dev_phi_raw[, 2] * sigma_phi_vec[2]
            dev_phi21 <- dev_phi_raw[, 3] * sigma_phi_vec[3]
            dev_phi22 <- dev_phi_raw[, 4] * sigma_phi_vec[4]

            # --- Generate Time Series Data for Each Subject ---
            subject_data_list <- vector("list", N_val)

            for (i_idx in seq_len(N_val)) {
                # Subject-specific parameters (will be global if RE SDs are 0)
                mu1_i <- mu_i[i_idx, 1]; mu2_i <- mu_i[i_idx, 2]
                phi11_i <- phi11_base + dev_phi11[i_idx]
                phi12_i <- phi12_base + dev_phi12[i_idx]
                phi21_i <- phi21_base + dev_phi21[i_idx]
                phi22_i <- phi22_base + dev_phi22[i_idx]

                y_mat <- matrix(NA_real_, nrow = T_total, ncol = 2)
                # Initialize first point at subject mean
                y_mat[1, ] <- c(mu1_i, mu2_i)

                # Generate time series
                for (t_idx in 2:T_total) {
                    y_prev <- y_mat[t_idx - 1, ]
                    if(anyNA(y_prev)) stop("NA encountered in y_prev, check initialization/loop logic.")

                    centered_prev <- c(y_prev[1] - mu1_i, y_prev[2] - mu2_i)
                    pred_1 <- phi11_i * centered_prev[1] + phi12_i * centered_prev[2]
                    pred_2 <- phi21_i * centered_prev[1] + phi22_i * centered_prev[2]

                    # Draw correlated innovations using copula
                    U <- draw_copula(n = 1, copula_info)
                    e1 <- draw_from_margin(U[1, 1], margin1_info)
                    e2 <- draw_from_margin(U[1, 2], margin2_info)

                    y_mat[t_idx, 1] <- mu1_i + pred_1 + e1
                    y_mat[t_idx, 2] <- mu2_i + pred_2 + e2
                } # End time loop

                # Discard burn-in period
                y_eff <- y_mat[(burnin + 1):T_total, , drop = FALSE]
                if (nrow(y_eff) != T_val) stop("Effective data length after burn-in is incorrect.")

                # Store subject data
                subject_data_list[[i_idx]] <- data.frame(
                    i = i_idx,
                    t = seq_len(T_val),
                    y1 = y_eff[, 1],
                    y2 = y_eff[, 2]
                )
            } # End subject loop

            # Combine data for all subjects for this replication
            big_df <- dplyr::bind_rows(subject_data_list)

            # --- Prepare Output List for this Replication ---
            true_params_list <- list(
                condition_id = current_cond_id,
                N = N_val, T = T_val,
                mu_global = mu_global_re_mean,
                phi_base = c(phi11_base, phi12_base, phi21_base, phi22_base),
                sigma_mu_vec = sigma_mu_vec, sigma_phi_vec = sigma_phi_vec,
                margin1 = margin1_info, margin2 = margin2_info,
                copula = copula_info,
                skew_normal_true_mean = 0.0, # True mean is zero by construction
                # Also save the factor levels for this condition for reference
                target_alpha = alpha_val, # Use saved alpha_val
                has_random_effects = cond$has_random_effects,
                data_seed = data_seed
            )

            sim_dat <- list(
                condition_id = current_cond_id,
                rep_i = rep_i,
                N = N_val, T = T_val,
                data = big_df,
                true_params = true_params_list
            )

            # --- Save Replication Data ---
            out_file <- file.path(output_dir, sprintf("sim_data_cond%03d_rep%03d.rds", current_cond_id, rep_i))
            saveRDS(sim_dat, file = out_file)

        } # End replication loop

        cat(sprintf("--- Finished Condition ID %d ---\n", current_cond_id))

    } # End condition loop

    cat("\nMultilevel VAR(1) simulation for all conditions completed.\n")
    invisible(NULL)
}