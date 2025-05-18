###########################################################################
# process_fits.R
#
# Intermediate script to process raw Stan fit objects.
# Reads fits from ../fits, simulation data/conditions from ../data.
# Extracts parameter summaries, sampler diagnostics, calculates metrics.
# Saves aggregated results required by evaluation.qmd and convergence.qmd
# into the current directory (assumed to be 'results').
#
# Run this script *after* fit_models.R has completed successfully
# and *before* running the Quarto reports.
###########################################################################

# --- Libraries ---
library(rstan)
library(tidyverse)
library(posterior) # Optional, but useful for some summaries
library(this.path) # To determine script location

# --- Configuration ---
cat("--- Starting Fit Processing ---\n")

# Set Base Directory relative to this script's location
RESULTS_DIR <- this.dir()
setwd(RESULTS_DIR) # Set WD to results directory for relative paths

# Define relative paths from the results directory
FITS_DIR <- file.path(RESULTS_DIR, "../fits")
DATA_DIR <- file.path(RESULTS_DIR, "../data")

# Check if directories exist
if (!dir.exists(FITS_DIR)) stop("Fits directory not found: ", normalizePath(FITS_DIR))
if (!dir.exists(DATA_DIR)) stop("Data directory not found: ", normalizePath(DATA_DIR))

# --- Load Simulation Conditions ---
sim_conds_file <- file.path(DATA_DIR, "sim_conditions.rds")
if (!file.exists(sim_conds_file)) stop("Cannot find sim conditions file: ", sim_conds_file)
sim_conditions_df <- readRDS(sim_conds_file) %>%
  # Ensure condition_id is integer for joining later
  mutate(condition_id = as.integer(condition_id))

cat("Loaded", nrow(sim_conditions_df), "simulation conditions.\n")

# --- Helper Functions ---

# Safely read RDS, checking if it's a stanfit object and has samples
safe_read_stanfit <- function(filename) {
  full_path <- filename
  if (!file.exists(full_path)) {
    warning(sprintf("File does not exist: %s", full_path), call. = FALSE)
    return(NULL)
  }
  fit <- NULL
  tryCatch({
    fit <- readRDS(full_path)
  }, error = function(e) {
    warning(sprintf("Error reading RDS file %s: %s", basename(full_path), conditionMessage(e)), call. = FALSE)
    return(NULL)
  })
  if (!is.null(fit)) {
    if (!inherits(fit, "stanfit")) {
      warning(sprintf("File %s is not a stanfit object.", basename(full_path)), call. = FALSE)
      return(NULL)
    }
    # Basic check for non-empty samples slot
    samples_exist <- FALSE
    tryCatch({
      if (!is.null(slot(fit, "sim")) && !is.null(slot(fit, "sim")$samples) &&
          length(slot(fit, "sim")$samples) > 0 && length(slot(fit, "sim")$samples[[1]]) > 0) {
        samples_exist <- TRUE
      }
    }, error = function(e_sample) { }) # Ignore errors accessing slot here
    
    if (!samples_exist) {
      warning(sprintf("Fit object seems empty or malformed (no samples): %s", basename(full_path)), call. = FALSE)
      return(NULL)
    }
    return(fit)
  } else {
    return(NULL)
  }
}


# E-FMI calculation (from evaluation.qmd)
calc_e_fmi <- function(energy_vec) {
  n <- length(energy_vec); if (n < 3) return(NA_real_)
  # Ensure finite values only
  energy_vec <- energy_vec[is.finite(energy_vec)]
  n <- length(energy_vec); if (n < 3) return(NA_real_)
  Emean <- mean(energy_vec); varE <- sum((energy_vec - Emean)^2) / (n - 1)
  # Calculate lag-1 autocorrelation robustly
  if (n < 2) return(NA_real_)
  lag1_sum <- sum((energy_vec[-1] - Emean) * (energy_vec[-n] - Emean))
  if (!is.finite(lag1_sum)) return(NA_real_)
  lag1 <- lag1_sum / (n - 1)
  denom <- varE + 2 * lag1; if (is.na(denom) || denom <= 1e-9) return(NA_real_) # Avoid division by zero/small
  efmi_val <- varE / denom
  # Clamp unreasonable values potentially caused by high noise
  pmax(0, pmin(1.5, efmi_val)) # Clamp between 0 and 1.5 (theoretical max > 1 possible but rare)
}


# Map Stan parameter names to true parameter values from the sim data list
get_true_value <- function(param_name, true_params_list, model_type) {
  true_val <- NA_real_
  # Handle indexed parameters first
  if (grepl("mu_global\\[\\d+\\]", param_name)) {
    idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", param_name))
    true_val <- true_params_list$mu_global[idx]
  } else if (grepl("sigma_mu\\[\\d+\\]", param_name)) {
    idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", param_name))
    true_val <- true_params_list$sigma_mu_vec[idx]
  } else if (grepl("sigma_phi\\[\\d+\\]", param_name)) {
    idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", param_name))
    true_val <- true_params_list$sigma_phi_vec[idx]
  } else if (grepl("xi\\[\\d+\\]", param_name)) {
    idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", param_name))
    margin_list <- if(idx == 1) true_params_list$margin1 else true_params_list$margin2
    true_val <- margin_list$params$xi
  } else if (grepl("omega\\[\\d+\\]", param_name)) {
    idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", param_name))
    margin_list <- if(idx == 1) true_params_list$margin1 else true_params_list$margin2
    true_val <- margin_list$params$omega
  } else if (grepl("alpha\\[\\d+\\]", param_name)) {
    idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", param_name))
    margin_list <- if(idx == 1) true_params_list$margin1 else true_params_list$margin2
    true_val <- margin_list$params$alpha
  } else if (grepl("sigma\\[\\d+\\]", param_name)) { # Normal model residual SD
    true_val <- NA_real_ # No direct equivalent in SN simulation (target variance was set)
    # Handle non-indexed parameters
  } else {
    true_val <- switch(param_name,
                       "phi11_base" = true_params_list$phi_base[1], # ML model name
                       "phi12_base" = true_params_list$phi_base[2], # ML model name
                       "phi21_base" = true_params_list$phi_base[3], # ML model name
                       "phi22_base" = true_params_list$phi_base[4], # ML model name
                       "phi11" = true_params_list$phi_base[1],      # SL model name maps to base
                       "phi12" = true_params_list$phi_base[2],      # SL model name maps to base
                       "phi21" = true_params_list$phi_base[3],      # SL model name maps to base
                       "phi22" = true_params_list$phi_base[4],      # SL model name maps to base
                       "rho" = true_params_list$copula$params$rho,
                       # Add cases for other specific parameters if needed
                       NA_real_ # Default if no match
    )
  }
  return(true_val)
}

# Assign parameter category based on name
get_param_category <- function(param_name) {
  if (grepl("^(mu_global|phi|rho)", param_name)) {
    "Fixed Effects / Copula"
  } else if (grepl("sigma_mu\\[\\d+\\]|sigma_phi\\[\\d+\\]", param_name)) {
    "RE SDs"
  } else if (grepl("var_sigma_mu\\[\\d+\\]|var_sigma_phi\\[\\d+\\]", param_name)) {
    "RE Variances"
  } else if (grepl("^(sigma|xi|omega|alpha)\\[\\d+\\]", param_name)) {
    "Residual Params"
  } else {
    "Other" # Catch-all
  }
}


# --- Find Fit Files ---
fit_files <- list.files(FITS_DIR, pattern = "^fit_.*_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)

if (length(fit_files) == 0) {
  stop("No fit files found in: ", normalizePath(FITS_DIR), ". Did the fitting step run and save files?")
}
cat("Found", length(fit_files), "potential fit files.\n")

# --- Process Fits in a Loop ---
param_results_list <- list()
sampler_results_list <- list()
files_processed <- 0
files_failed <- 0

# Use purrr::map or a loop
for (fit_file in fit_files) {
  cat(".") # Progress indicator
  base_name <- basename(fit_file)
  
  # Extract info from filename using regex
  match_info <- stringr::str_match(base_name, "^fit_(normal|skewnormal)_cond(\\d+)_rep(\\d+)\\.rds$")
  if (is.na(match_info[1,1])) {
    warning("Could not parse filename: ", base_name, call. = FALSE)
    files_failed <- files_failed + 1
    next
  }
  model_type <- match_info[1, 2]
  cond_id <- as.integer(match_info[1, 3])
  rep_i <- as.integer(match_info[1, 4])
  
  # --- Load Stan Fit Object ---
  fit <- safe_read_stanfit(fit_file)
  if (is.null(fit)) {
    warning("Skipping invalid or unreadable fit file: ", base_name, call. = FALSE)
    files_failed <- files_failed + 1
    next
  }
  
  # --- Load Corresponding Simulation Data to get True Parameters ---
  sim_data_file <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cond_id, rep_i))
  if (!file.exists(sim_data_file)) {
    warning("Simulation data file not found for ", base_name, ": ", sim_data_file, call. = FALSE)
    warning("Cannot determine true parameter values for this fit. Skipping.", call. = FALSE)
    files_failed <- files_failed + 1
    next
  }
  sim_data <- readRDS(sim_data_file)
  true_params <- sim_data$true_params
  if (is.null(true_params)) {
    warning("True parameters missing in simulation data file: ", basename(sim_data_file), call. = FALSE)
    warning("Cannot determine true parameter values for this fit. Skipping.", call. = FALSE)
    files_failed <- files_failed + 1
    next
  }
  
  
  # --- Extract Parameter Summaries ---
  summary_df <- tryCatch({
    summary(fit)$summary %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      as_tibble() %>%
      # Select and rename relevant columns
      select(parameter,
             post_mean = mean, post_median = `50%`, post_sd = sd,
             ci_low = `2.5%`, ci_high = `97.5%`,
             n_eff, Rhat)
  }, error = function(e) {
    warning("Error getting summary for ", base_name, ": ", conditionMessage(e), call. = FALSE)
    return(NULL)
  })
  
  if (is.null(summary_df)) {
    files_failed <- files_failed + 1
    next # Skip to next file if summary failed
  }
  
  # --- Add True Values, Bias, Coverage ---
  summary_df <- summary_df %>%
    mutate(
      condition_id = cond_id,
      rep_i = rep_i,
      model_type = model_type,
      # Use rowwise apply or map to get true value for each parameter
      true_value = map_dbl(parameter, ~get_true_value(.x, true_params, model_type)),
      bias = post_mean - true_value,
      coverage = ifelse(!is.na(true_value), (true_value >= ci_low & true_value <= ci_high), NA),
      ci_width = ci_high - ci_low,
      rel_bias = ifelse(abs(true_value) > 1e-6 & !is.na(bias), bias / abs(true_value), NA_real_) # Avoid division by zero
    )
  
  # --- Add Derived Variance Parameters ---
  # Check if sigma_mu/sigma_phi exist before trying to square them
  if(all(paste0("sigma_mu[", 1:2, "]") %in% summary_df$parameter)) {
    var_mu_df <- summary_df %>% filter(grepl("sigma_mu\\[\\d\\]", parameter)) %>%
      mutate(var_estimate = post_mean^2, # Simple estimate of variance
             parameter = gsub("sigma_mu", "var_sigma_mu", parameter),
             true_value = true_value^2, # Transform true value too
             # Recalculate metrics based on variance
             bias = var_estimate - true_value,
             # Note: Coverage/CI width for variance requires recalculation from posterior draws,
             #       this is just a point estimate transformation. Set related metrics to NA.
             post_mean = var_estimate, post_median=NA, post_sd=NA, ci_low=NA, ci_high=NA, coverage=NA, ci_width=NA, rel_bias = NA
      ) %>% select(-var_estimate)
    summary_df <- bind_rows(summary_df, var_mu_df)
  }
  if(all(paste0("sigma_phi[", 1:4, "]") %in% summary_df$parameter)) {
    var_phi_df <- summary_df %>% filter(grepl("sigma_phi\\[\\d\\]", parameter)) %>%
      mutate(var_estimate = post_mean^2,
             parameter = gsub("sigma_phi", "var_sigma_phi", parameter),
             true_value = true_value^2,
             bias = var_estimate - true_value,
             post_mean = var_estimate, post_median=NA, post_sd=NA, ci_low=NA, ci_high=NA, coverage=NA, ci_width=NA, rel_bias = NA
      ) %>% select(-var_estimate)
    summary_df <- bind_rows(summary_df, var_phi_df)
  }
  
  # --- Add Parameter Category ---
  summary_df <- summary_df %>%
    mutate(param_category = map_chr(parameter, get_param_category))
  
  param_results_list[[length(param_results_list) + 1]] <- summary_df
  
  
  # --- Extract Sampler Diagnostics ---
  sampler_params <- tryCatch({
    rstan::get_sampler_params(fit, inc_warmup = FALSE)
  }, error = function(e) {
    warning("Error getting sampler params for ", base_name, ": ", conditionMessage(e), call. = FALSE)
    return(NULL)
  })
  
  if (!is.null(sampler_params)) {
    # Summarize across chains and iterations
    # Convert list of matrices to a tidy data frame
    sampler_params_df <- map_dfr(seq_along(sampler_params), ~as_tibble(sampler_params[[.x]]) %>% mutate(chain = .x))
    
    # Calculate summary stats per fit
    total_divergences <- sum(sampler_params_df$divergent__, na.rm = TRUE)
    maxdepth_exceeded <- sum(sampler_params_df$treedepth__ >= attr(fit@sim$samples[[1]], "sampler_params")$max_treedepth, na.rm = TRUE) # Max treedepth might vary
    
    # Calculate E-FMI (average across chains)
    eFMI_per_chain <- map_dbl(sampler_params, ~calc_e_fmi(.x[, "energy__"]))
    avg_eFMI <- mean(eFMI_per_chain, na.rm = TRUE) # Average E-FMI across chains
    
    sampler_summary <- tibble(
      condition_id = cond_id,
      rep_i = rep_i,
      model_type = model_type,
      divergences = total_divergences,
      maxdepth_exceeded = maxdepth_exceeded,
      eFMI = avg_eFMI # Store average E-FMI
    )
    sampler_results_list[[length(sampler_results_list) + 1]] <- sampler_summary
  } else {
    # Add placeholder if sampler params failed
    sampler_results_list[[length(sampler_results_list) + 1]] <- tibble(
      condition_id = cond_id, rep_i = rep_i, model_type = model_type,
      divergences = NA_integer_, maxdepth_exceeded = NA_integer_, eFMI = NA_real_
    )
  }
  
  files_processed <- files_processed + 1
  
} # End loop over fit files
cat("\n") # Newline after progress dots

cat("Processed", files_processed, "fit files.\n")
if (files_failed > 0) {
  warning("Failed to process ", files_failed, " files (see warnings above).", call. = FALSE)
}


# --- Aggregate Results ---
if (length(param_results_list) > 0) {
  all_param_results <- bind_rows(param_results_list)
  
  # --- Merge with Simulation Condition Factors ---
  # Select only the factor columns from sim_conditions_df to avoid duplicating metrics
  sim_factors <- sim_conditions_df %>%
    select(condition_id, N, T, target_alpha, phi11_base, phi12_base, has_random_effects)
  
  all_param_results_merged <- all_param_results %>%
    left_join(sim_factors, by = "condition_id")
  
  # --- Save Parameter Results ---
  param_output_file <- file.path(RESULTS_DIR, "simulation_parameter_results_with_vars.rds")
  saveRDS(all_param_results_merged, file = param_output_file)
  cat("Saved parameter results to:", param_output_file, "\n")
  cat("Total parameters summarised:", nrow(all_param_results_merged), "\n")
  
} else {
  warning("No parameter results were successfully processed.", call. = FALSE)
  all_param_results_merged <- data.frame() # Ensure it exists
}

if (length(sampler_results_list) > 0) {
  all_sampler_results <- bind_rows(sampler_results_list)
  
  # --- Merge with Simulation Condition Factors ---
  # Re-use sim_factors defined above
  all_sampler_results_merged <- all_sampler_results %>%
    left_join(sim_factors, by = "condition_id")
  
  # --- Save Sampler Info ---
  sampler_output_file <- file.path(RESULTS_DIR, "sampler_information.rds")
  saveRDS(all_sampler_results_merged, file = sampler_output_file)
  cat("Saved sampler information to:", sampler_output_file, "\n")
  cat("Total fits summarised for sampler info:", nrow(all_sampler_results_merged), "\n")
} else {
  warning("No sampler results were successfully processed.", call. = FALSE)
  all_sampler_results_merged <- data.frame() # Ensure it exists
}


cat("--- Fit Processing Finished ---\n")