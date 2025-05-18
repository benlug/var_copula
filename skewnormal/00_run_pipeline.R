###########################################################################
# run_pipeline.R
#
# Main control script for the multilevel VAR(1) Skew-Normal simulation study.
# Defines simulation conditions, runs simulation, checks, fitting, and analysis.
###########################################################################

# --- Libraries ---
# tidyverse helpers
library(dplyr)
library(tidyr)
library(purrr)
# to locate the script itself
library(this.path)
# compile / run stan file
library(rstan)

# --- Configuration ---
BASE_DIR <- this.dir() # this resolves the folder that contains this file
setwd(BASE_DIR)

# --- Simulation Factors & Levels ---
# factor_levels <- list(
#     skew_alpha = c(-10.0, 10.0),  # Negative (left) and Positive (right) skew alpha
#     T_levels = c(30, 100),        # Time points per subject
#     N_levels = c(50, 300),        # Number of subjects
#     ar_coeff = c(0.2, 0.5),       # Base diagonal AR(1) coefficient (phi11=phi22)
#     cross_lag_coeff = c(0.0, 0.3), # Base cross-lag coefficient (phi12=phi21)
#     random_effects = c(TRUE, FALSE), # Include random intercepts and slopes?
#     replications = 1 #
# )
factor_levels <- list(
    skew_alpha = c(-10.0, 10.0), # direction and magnitude of marginal skew
    T_levels = c(30, 100), # time points per subject
    N_levels = c(50), # number of subjects / persons
    ar_coeff = c(0.2), # ar coeffs ϕ11 = ϕ22
    cross_lag_coeff = c(0.1), # cross-lag ϕ12 = ϕ21
    random_effects = c(TRUE), # multilevel VAR (random µ and ϕ)
    replications = 8 # monte carlo reps per condition
)

# --- Fixed Simulation Parameters ---
# All subjects share these unless random_effects = FALSE, in which case the SDs are replaced with zeros later on
fixed_params <- list(
    target_variance = 1, # marginal variance of each univariate margin
    copula_rho = 0.4, # gaussian copula residual correlation
    mu_global_mean = 0.0, # grand mean of random intercepts
    # Standard deviations for random effects IF random_effects = TRUE
    # Set to 0 if random_effects = FALSE within the scenario generation
    sigma_mu_sd_if_re = c(0.5, 0.5), # SD of random intercepts (µ)
    sigma_phi_sd_if_re = c(0.1, 0.1, 0.1, 0.1), # SD of random VAR coeffs (ϕ’s)
    burnin = 50 # discard these points before saving data
)

# --- Pipeline Control Flags ---
# Set these to TRUE/FALSE to run specific parts of the pipeline
run_simulation <- TRUE # generate raw data
run_checks <- TRUE # runs visual checks on simulated data
run_fitting <- TRUE # fit stan models
run_processing <- TRUE # extract estimates into tidy tables
run_evaluation <- TRUE # render Quarto parameter-recovery report
run_convergence <- TRUE # render Quarto convergence report

# --- Directory Setup ---
DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULTS_DIR <- file.path(BASE_DIR, "results")

if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR)
if (!dir.exists(CHECKS_DIR)) dir.create(CHECKS_DIR)
if (!dir.exists(FITS_DIR)) dir.create(FITS_DIR)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR)

# --- Generate Simulation Conditions ---
var_scenarios_grid <- expand.grid(
    target_alpha = factor_levels$skew_alpha,
    T = factor_levels$T_levels,
    N = factor_levels$N_levels,
    phi11_base = factor_levels$ar_coeff,
    phi12_base = factor_levels$cross_lag_coeff,
    has_random_effects = factor_levels$random_effects
)
var_scenarios_grid$phi22_base <- var_scenarios_grid$phi11_base
var_scenarios_grid$phi21_base <- var_scenarios_grid$phi12_base

# add fixed parameters and calculate derived ones
var_scenarios <- var_scenarios_grid %>%
    mutate(
        n_reps = factor_levels$replications,
        mu_global_1 = fixed_params$mu_global_mean,
        mu_global_2 = fixed_params$mu_global_mean,
        target_variance = fixed_params$target_variance,
        copula_rho = fixed_params$copula_rho,
        burnin = fixed_params$burnin,
        # set RE SDs based on has_random_effects flag
        sigma_mu_vec = map(has_random_effects, ~ if (.x) fixed_params$sigma_mu_sd_if_re else c(0.0, 0.0)),
        sigma_phi_vec = map(has_random_effects, ~ if (.x) fixed_params$sigma_phi_sd_if_re else c(0.0, 0.0, 0.0, 0.0)),
        # --- Calculate Skew-Normal params (xi, omega) for zero mean ---
        # this calculation depends only on alpha and target variance
        delta_val = target_alpha / sqrt(1 + target_alpha^2),
        denom_var = (1 - (2 * delta_val^2 / pi)),
        # add check for invalid parameters
        required_omega_sq = ifelse(denom_var > 1e-9, target_variance / denom_var, NA_real_),
        target_omega = ifelse(!is.na(required_omega_sq) & required_omega_sq > 0, sqrt(required_omega_sq), NA_real_),
        required_xi = ifelse(!is.na(target_omega), -target_omega * delta_val * sqrt(2 / pi), NA_real_)
    ) %>%
    filter(!is.na(target_omega)) # Remove conditions where variance isn't achievable

# check if any invalid scenarios were created
if (nrow(var_scenarios) == 0) {
    stop("No valid simulation scenarios generated. Check fixed_params$target_variance and factor_levels$skew_alpha.")
}

# create margin and copula lists required by simulate_data.R
var_scenarios <- var_scenarios %>%
    mutate(
        margin1 = pmap(
            list(xi = required_xi, omega = target_omega, alpha = target_alpha),
            ~ list(dist = "skewnormal", params = list(xi = ..1, omega = ..2, alpha = ..3))
        ),
        margin2 = pmap(
            list(xi = required_xi, omega = target_omega, alpha = target_alpha),
            ~ list(dist = "skewnormal", params = list(xi = ..1, omega = ..2, alpha = ..3))
        ),
        copula_choice = map(copula_rho, ~ list(copula = "gaussian", params = list(rho = .x)))
    )

# Select final columns for the simulation conditions file
# Note: Added has_random_effects and original factor columns for clarity in results
sim_conditions_to_save <- var_scenarios %>%
    select(
        N, T, mu_global_1, mu_global_2,
        phi11_base, phi12_base, phi21_base, phi22_base,
        sigma_mu_vec, sigma_phi_vec, n_reps, burnin,
        margin1, margin2, copula_choice,
        # Keep original factors for easier analysis later
        target_alpha, has_random_effects
    ) %>%
    # Add a unique condition index
    mutate(condition_id = row_number()) %>%
    select(condition_id, everything())


# Save the conditions
conditions_file <- file.path(DATA_DIR, "sim_conditions.rds")
saveRDS(sim_conditions_to_save, conditions_file)
cat(sprintf("Generated %d simulation conditions.\n", nrow(sim_conditions_to_save)))
cat("Saved simulation conditions to:", conditions_file, "\n")

# --- Execute Pipeline Steps ---

# 1. Simulation
if (run_simulation) {
    cat("\n--- Running Data Simulation ---\n")
    # Source the simulation script, which now contains the main function
    source(file.path(BASE_DIR, "simulate_data.R"), local = TRUE) # local=TRUE to keep functions local
    # Call the main simulation function from the sourced script
    simulate_all_conditions(
        sim_conditions_df = sim_conditions_to_save,
        output_dir = DATA_DIR
    )
    cat("--- Data Simulation Finished ---\n")
} else {
    cat("\n--- Skipping Data Simulation ---\n")
}

# 2. Simulation Checks
if (run_checks) {
    cat("\n--- Running Simulation Checks ---\n")
    source(file.path(BASE_DIR, "check_simulations.R"), local = TRUE)
    # The check script should read conditions and loop through data files
    run_post_sim_checks(
        data_dir = DATA_DIR,
        checks_dir = CHECKS_DIR
    )
    cat("--- Simulation Checks Finished ---\n")
} else {
    cat("\n--- Skipping Simulation Checks ---\n")
}

# 3. Model Fitting
if (run_fitting) {
    cat("\n--- Running Model Fitting ---\n")
    # Source the modified fit_models.R which now handles parallelization
    source(file.path(BASE_DIR, "fit_models.R"), local = TRUE)
    # Call the main fitting function (it will run in parallel internally)
    if (TRUE) {
        fit_multilevel_models_sn_vec(
            data_dir = DATA_DIR,
            fits_dir = FITS_DIR,
            stan_models_dir = BASE_DIR
        ) # Pass base dir for .stan files
        # Optional: add num_cores argument here if desired
        # num_cores = 4)
    }

    # DEBUGGING CALL:
    # if (FALSE) {
    #     fit_multilevel_models_sn_vec(
    #         data_dir = DATA_DIR,
    #         fits_dir = FITS_DIR,
    #         stan_models_dir = BASE_DIR,
    #         debug_mode = TRUE, # KEEP THIS TRUE
    #         debug_file_limit = 5, # Keep this set
    #         stan_chains = 4 # Ensure this matches cores you want for chains
    #     )
    # }
    cat("--- Model Fitting Finished ---\n")
} else {
    cat("\n--- Skipping Model Fitting ---\n")
}


# 3.5 Process Fits
if (run_processing) {
    cat("\n--- Running Post-Fit Processing ---\n")
    process_script <- file.path(RESULTS_DIR, "process_fits.R")
    if (file.exists(process_script)) {
        tryCatch(
            {
                source(process_script, local = TRUE) # Source locally
                cat("--- Post-Fit Processing Finished ---\n")
            },
            error = function(e) {
                cat("ERROR running process_fits.R:", conditionMessage(e), "\n")
                stop("Post-processing failed.") # Stop pipeline if processing fails
            }
        )
    } else {
        cat("WARNING: process_fits.R not found in results directory. Cannot generate summary files.\n")
        stop("Processing script missing.") # Stop pipeline if script is missing
    }
} else {
    cat("\n--- Skipping Post-Fit Processing ---\n")
    # Optional: Check if files exist anyway if skipping processing?
    # Might be needed if reports are run independently later.
    param_results_file_check <- file.path(RESULTS_DIR, "simulation_parameter_results_with_vars.rds")
    sampler_info_file_check <- file.path(RESULTS_DIR, "sampler_information.rds")
    if (!file.exists(param_results_file_check) || !file.exists(sampler_info_file_check)) {
        warning("Skipping processing, but required summary files are missing. Reports may fail.", call. = FALSE)
    }
}


# --- Analysis Steps (Quarto Reports) ---
# Ensure output directories for reports are set correctly *inside* the qmd files
# or pass them as execution parameters if needed. Here, we assume they output
# relative to their own location or to RESULTS_DIR.

# 4. Parameter Recovery Evaluation
if (run_evaluation) {
    cat("\n--- Running Parameter Recovery Evaluation ---\n")
    eval_qmd_file <- file.path(BASE_DIR, "evaluation.qmd")
    if (file.exists(eval_qmd_file)) {
        tryCatch(
            {
                quarto_render(
                    input = eval_qmd_file,
                    # Execute parameters can pass paths to the qmd file
                    # execute_params = list(fits_dir = FITS_DIR,
                    #                       data_dir = DATA_DIR,
                    #                       results_dir = RESULTS_DIR),
                    output_file = file.path(RESULTS_DIR, "parameter_evaluation_report.html")
                )
                cat("Rendered evaluation report.\n")
            },
            error = function(e) {
                cat("ERROR rendering evaluation.qmd:", conditionMessage(e), "\n")
            }
        )
    } else {
        cat("WARNING: evaluation.qmd not found at:", eval_qmd_file, "\n")
    }
    cat("--- Parameter Recovery Evaluation Finished ---\n")
} else {
    cat("\n--- Skipping Parameter Recovery Evaluation ---\n")
}

# 5. Convergence Diagnostics
if (run_convergence) {
    cat("\n--- Running Convergence Diagnostics ---\n")
    conv_qmd_file <- file.path(BASE_DIR, "convergence.qmd")
    if (file.exists(conv_qmd_file)) {
        tryCatch(
            {
                quarto_render(
                    input = conv_qmd_file,
                    # execute_params = list(fits_dir = FITS_DIR,
                    #                       data_dir = DATA_DIR,
                    #                       results_dir = RESULTS_DIR),
                    output_file = file.path(RESULTS_DIR, "convergence_diagnostics_report.html")
                )
                cat("Rendered convergence report.\n")
            },
            error = function(e) {
                cat("ERROR rendering convergence.qmd:", conditionMessage(e), "\n")
            }
        )
    } else {
        cat("WARNING: convergence.qmd not found at:", conv_qmd_file, "\n")
    }
    cat("--- Convergence Diagnostics Finished ---\n")
} else {
    cat("\n--- Skipping Convergence Diagnostics ---\n")
}

cat("\n--- Pipeline execution complete. ---\n")
