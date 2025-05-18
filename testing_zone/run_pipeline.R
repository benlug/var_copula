###########################################################################
# run_pipeline.R
#
# Main control script for the single-level VAR(1) Copula simulation study.
# Allows for DIFFERENT alpha values for the two marginal distributions.
# Simplified naming and reduced verbosity.
###########################################################################

# --- Libraries ---
library(dplyr)
library(tidyr)
library(purrr)
library(this.path)
library(rstan)
library(copula)

# --- Configuration ---
BASE_DIR <- this.dir()
setwd(BASE_DIR)

# --- Simulation Factors & Levels ---
# factor_levels <- list(
#   dgp_copula = c("gaussian", "clayton"),
#   dgp_alpha1 = c(-5.0, 0.0, 5.0), # Alpha for margin 1
#   dgp_alpha2 = c(-5.0, 0.0, 5.0), # Alpha for margin 2
#   dgp_tau = c(0.2, 0.5),
#   T_levels = c(30, 100),
#   phi11 = c(0.2, 0.6),
#   phi22 = c(0.2, 0.6),
#   phi12 = c(0.0, 0.2),
#   phi21 = c(0.0, 0.2),
#   replications = 100 # SET TO SMALL NUMBER FOR TESTING (e.g., 2-4)
#   # replications = 4
# )

# --- Simulation Factors & Levels ---
# factor_levels <- list(
#   dgp_copula = c("clayton"),
#   dgp_alpha1 = c(-9), # Alpha for margin 1
#   dgp_alpha2 = c(9), # Alpha for margin 2
#   dgp_tau = c(0.5),
#   T_levels = c(30, 100),
#   phi11 = c(0.5),
#   phi22 = c(0.5),
#   phi12 = c(0.2),
#   phi21 = c(0.2),
#   replications = 32 # SET TO SMALL NUMBER FOR TESTING (e.g., 2-4)
#   # replications = 4
# )

factor_levels <- list(
  dgp_copula = c("gaussian", "clayton"),
  dgp_alpha1 = c(-5.0, 5.0),
  dgp_alpha2 = c(-5.0, 5.0),
  dgp_tau = c(0.2, 0.5),
  T_levels = c(30, 100),
  phi11 = c(0.6),
  phi22 = c(0.6),
  phi12 = c(0.2),
  phi21 = c(0.2),
  replications = 24
)


# --- Fixed Simulation Parameters ---
fixed_params <- list(
  target_variance = 1.0,
  mu_intercepts = c(0.0, 0.0),
  burnin = 30
)

# --- Pipeline Control Flags ---
run_simulation <- TRUE
run_checks <- TRUE
run_fitting <- TRUE
run_processing <- FALSE # Assumes process_fits.R is updated for separate alphas
# run_evaluation <- FALSE # Requires updated Quarto/R script
# run_convergence <- FALSE # Requires updated Quarto/R script

# --- Directory Setup (Simplified Names) ---
DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULTS_DIR <- file.path(BASE_DIR, "results")

# Create directories if they don't exist
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR)
if (!dir.exists(CHECKS_DIR)) dir.create(CHECKS_DIR)
if (!dir.exists(FITS_DIR)) dir.create(FITS_DIR)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR)

# --- Helper Functions ---
tau_to_copula_param <- function(tau, copula_family) {
  if (copula_family == "gaussian") {
    cop <- normalCopula(param = 0.5, dim = 2)
    param <- iTau(cop, tau)
    param <- max(-0.9999, min(0.9999, param))
    return(list(type = "gaussian", param_name = "rho", value = param))
  } else if (copula_family == "clayton") {
    cop <- claytonCopula(param = 0.5, dim = 2)
    param <- iTau(cop, tau)
    param <- max(1e-6, param)
    return(list(type = "clayton", param_name = "theta", value = param))
  } else {
    stop("Unsupported copula family: ", copula_family)
  }
}

calculate_sn_params <- function(alpha, target_variance = 1.0) {
  is_normal <- (abs(alpha) < 1e-6)
  if (is_normal) {
    return(list(xi = 0, omega = sqrt(target_variance), alpha = 0))
  } else {
    delta_val <- alpha / sqrt(1 + alpha^2)
    denom_var <- (1 - (2 * delta_val^2 / pi))
    required_omega_sq <- ifelse(denom_var > 1e-9, target_variance / denom_var, NA_real_)
    target_omega <- ifelse(!is.na(required_omega_sq) & required_omega_sq > 0, sqrt(required_omega_sq), NA_real_)
    required_xi <- ifelse(!is.na(target_omega), -target_omega * delta_val * sqrt(2 / pi), NA_real_)
    if (is.na(target_omega) || is.na(required_xi)) stop("Could not calculate valid xi/omega for alpha=", alpha)
    return(list(xi = required_xi, omega = target_omega, alpha = alpha))
  }
}

# --- Generate Simulation Conditions ---
cat("Generating simulation condition grid...\n")
var_scenarios_grid <- expand.grid(
  dgp_copula_type = factor_levels$dgp_copula,
  dgp_alpha1 = factor_levels$dgp_alpha1,
  dgp_alpha2 = factor_levels$dgp_alpha2,
  dgp_tau = factor_levels$dgp_tau,
  T = factor_levels$T_levels,
  phi11 = factor_levels$phi11,
  phi12 = factor_levels$phi12,
  phi21 = factor_levels$phi21,
  phi22 = factor_levels$phi22,
  stringsAsFactors = FALSE
)

var_scenarios <- var_scenarios_grid %>%
  mutate(
    n_reps = factor_levels$replications,
    mu1 = fixed_params$mu_intercepts[1],
    mu2 = fixed_params$mu_intercepts[2],
    target_variance = fixed_params$target_variance,
    burnin = fixed_params$burnin
  ) %>%
  rowwise() %>%
  mutate(
    sn_params1 = list(calculate_sn_params(dgp_alpha1, target_variance)),
    dgp_margin1_dist = ifelse(abs(dgp_alpha1) < 1e-6, "normal", "skewnormal"),
    dgp_margin1_params = ifelse(dgp_margin1_dist == "normal",
      list(list(mean = 0, sd = sn_params1$omega)),
      list(list(xi = sn_params1$xi, omega = sn_params1$omega, alpha = sn_params1$alpha))
    ),
    sn_params2 = list(calculate_sn_params(dgp_alpha2, target_variance)),
    dgp_margin2_dist = ifelse(abs(dgp_alpha2) < 1e-6, "normal", "skewnormal"),
    dgp_margin2_params = ifelse(dgp_margin2_dist == "normal",
      list(list(mean = 0, sd = sn_params2$omega)),
      list(list(xi = sn_params2$xi, omega = sn_params2$omega, alpha = sn_params2$alpha))
    ),
    dgp_copula_info = list(tau_to_copula_param(dgp_tau, dgp_copula_type))
  ) %>%
  ungroup() %>%
  mutate(
    correct_fit_margin_type = ifelse(abs(dgp_alpha1) < 1e-6 & abs(dgp_alpha2) < 1e-6, "normal", "skewnormal"),
    correct_fit_copula_type = dgp_copula_type,
    correct_fit_model_code = paste0(
      toupper(substr(correct_fit_margin_type, 1, 1)),
      toupper(substr(correct_fit_copula_type, 1, 1))
    )
  ) %>%
  select(-sn_params1, -sn_params2)

if (nrow(var_scenarios) == 0) stop("No valid simulation scenarios generated.")

sim_conditions_to_save <- var_scenarios %>%
  select(
    dgp_copula_type, dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11, phi12, phi21, phi22,
    n_reps, mu1, mu2, burnin,
    dgp_margin1_dist, dgp_margin1_params,
    dgp_margin2_dist, dgp_margin2_params,
    dgp_copula_info, correct_fit_model_code
  ) %>%
  mutate(condition_id = row_number()) %>%
  select(condition_id, everything())

conditions_file <- file.path(DATA_DIR, "sim_conditions.rds") # Simplified name
saveRDS(sim_conditions_to_save, conditions_file)
cat(sprintf("Generated %d conditions. Saved to: %s\n", nrow(sim_conditions_to_save), conditions_file))

# --- Execute Pipeline Steps ---

# 1. Simulation
if (run_simulation) {
  cat("\n--- Running Simulation ---\n")
  source(file.path(BASE_DIR, "simulate_data.R"), local = TRUE) # Use simplified simulation script name
  simulate_all_conditions_var1( # Function name from simulate_data.R
    sim_conditions_df = sim_conditions_to_save,
    output_dir = DATA_DIR
  )
  cat("--- Simulation Finished ---\n")
} else {
  cat("\n--- Skipping Simulation ---\n")
}

# 2. Simulation Checks
if (run_checks) {
  cat("\n--- Running Simulation Checks ---\n")
  check_script_path <- file.path(BASE_DIR, "check_simulations.R") # Use simplified name
  if (!file.exists(check_script_path)) {
    warning("Check script not found: ", check_script_path, ". Skipping checks.", call. = FALSE)
  } else {
    source(check_script_path, local = TRUE)
    run_post_sim_checks_var1( # Function name from check_simulations.R
      data_dir = DATA_DIR,
      checks_dir = CHECKS_DIR
    )
    cat("--- Simulation Checks Finished ---\n")
  }
} else {
  cat("\n--- Skipping Simulation Checks ---\n")
}

# 3. Model Fitting
if (run_fitting) {
  cat("\n--- Running Model Fitting ---\n")
  fit_script_path <- file.path(BASE_DIR, "fit_models.R") # Use simplified name
  if (!file.exists(fit_script_path)) {
    stop("Fitting script not found: ", fit_script_path)
  } else {
    source(fit_script_path, local = TRUE)
    fit_var1_copula_models( # Function name from fit_models.R
      data_dir = DATA_DIR,
      fits_dir = FITS_DIR,
      stan_models_dir = BASE_DIR, # Assumes .stan files are in BASE_DIR
      sim_conditions_file = conditions_file,
      num_cores = parallel::detectCores() - 4 # Optional parallel
    )
    cat("--- Model Fitting Finished ---\n")
  }
} else {
  cat("\n--- Skipping Model Fitting ---\n")
}

# 3.5 Process Fits
if (run_processing) {
  cat("\n--- Running Post-Fit Processing ---\n")
  # Assumes process_fits.R exists in BASE_DIR and is updated
  process_script_path <- file.path(RESULTS_DIR, "process_fits.R")
  if (!file.exists(process_script_path)) {
    warning("Processing script 'process_fits.R' not found in RESULTS_DIR", call. = FALSE)
    stop("Processing script missing.")
  }
  warning("Executing process_fits.R - Ensure it handles separate alpha1/alpha2!", call. = FALSE)
  tryCatch(
    {
      source(process_script_path, local = TRUE)
      cat("--- Post-Fit Processing Finished ---\n")
    },
    error = function(e) {
      cat("ERROR running process_fits.R:", conditionMessage(e), "\n")
      stop("Post-processing failed.")
    }
  )
} else {
  cat("\n--- Skipping Post-Fit Processing ---\n")
  param_results_file_check <- file.path(RESULTS_DIR, "parameter_summary.rds") # Simplified name
  sampler_info_file_check <- file.path(RESULTS_DIR, "sampler_summary.rds") # Simplified name
  if (!file.exists(param_results_file_check) || !file.exists(sampler_info_file_check)) {
    warning("Skipping processing, but required summary files are missing.", call. = FALSE)
  }
}

# --- Analysis Steps (Require Updates) ---
cat("\n--- Analysis/Reporting steps skipped (require updated scripts) ---\n")

cat("\n--- Pipeline execution complete (Sim, Checks, Fit, Process). ---\n")
cat("--- NOTE: Analysis/Reporting scripts require updates. ---\n")
