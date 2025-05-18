##########################################################
# fit_models.R (Generalized - Enhanced Debugging)
#
# Fits multilevel OR single-level VAR(1) models using vectorized Stan code.
# Handles parallel processing of simulation files internally.
# Selects appropriate Stan model (ML vs SL) based on condition info.
# Includes debug mode for sequential execution and verbose output.
# Called by run_pipeline.R.
##########################################################

# --- Libraries ---
# Load necessary libraries here
library(rstan)
library(dplyr)
library(tidyr)
library(foreach) # Added for parallelization
library(doParallel) # Added for parallelization

# Null default operator (useful for parameters)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Helper Function for Fitting WORKER (with Debugging) ---
# This function performs the work for a SINGLE data file on a parallel worker or sequentially.
# It receives the correct compiled models (ML or SL) based on the condition.
fit_models_for_single_file_worker <- function(data_filepath,
                                              compiled_norm_model, # Will be the correct one (ML or SL)
                                              compiled_skew_model, # Will be the correct one (ML or SL)
                                              fits_dir,
                                              is_multilevel, # Flag indicating model type
                                              # Stan control params
                                              prior_sigma_mu_scale, # Needed only if is_multilevel=TRUE
                                              prior_sigma_phi_scale, # Needed only if is_multilevel=TRUE
                                              stan_iter, stan_warmup, stan_chains,
                                              stan_adapt_delta, stan_max_treedepth,
                                              debug_mode = FALSE # Added debug flag
) {
  # --- Packages needed on worker (if parallel) ---
  # Ensure these are loaded via .packages in foreach if running in parallel
  require(rstan)
  require(dplyr)
  require(tidyr)
  require(stats) # Explicitly require stats for rnorm, sd etc.

  # --- Parse filename ---
  fname_base <- basename(data_filepath)
  cond_match <- regexpr("(?<=_cond)(\\d+)", fname_base, perl = TRUE)
  rep_match <- regexpr("(?<=_rep)(\\d+)(?=\\.rds$)", fname_base, perl = TRUE)
  if (cond_match == -1 || attr(cond_match, "match.length") == -1 ||
    rep_match == -1 || attr(rep_match, "match.length") == -1) {
    return(paste("ERROR: Could not parse cond/rep from", fname_base))
  }
  cond_num <- tryCatch(as.integer(regmatches(fname_base, cond_match)), error = function(e) NA)
  rep_num <- tryCatch(as.integer(regmatches(fname_base, rep_match)), error = function(e) NA)
  if (is.na(cond_num) || is.na(rep_num)) {
    return(paste("ERROR: Failed integer conversion for cond/rep from", fname_base))
  }


  # --- Define Output Filenames ---
  out_file_sn <- file.path(fits_dir, sprintf("fit_skewnormal_cond%03d_rep%03d.rds", cond_num, rep_num))
  out_file_norm <- file.path(fits_dir, sprintf("fit_normal_cond%03d_rep%03d.rds", cond_num, rep_num))

  # --- Skip if BOTH fits already exist ---
  if (file.exists(out_file_sn) && file.exists(out_file_norm)) {
    return(paste("Skipped:", fname_base, "(both fits already exist)"))
  }

  # --- Start Debug Output ---
  if (debug_mode) {
    cat("\n=====================================================\n")
    cat(sprintf("DEBUG MODE: Processing file: %s\n", fname_base))
    cat(sprintf("Condition: %d, Replication: %d\n", cond_num, rep_num))
    cat(sprintf("Is Multilevel: %s\n", is_multilevel))
    cat(sprintf("Using Normal Model: %s\n", compiled_norm_model@model_name))
    cat(sprintf("Using SkewNormal Model: %s\n", compiled_skew_model@model_name))
    cat(sprintf("Output SN file: %s\n", out_file_sn))
    cat(sprintf("Output Norm file: %s\n", out_file_norm))
    cat("-----------------------------------------------------\n")
  }

  # --- Load and Prepare Data ---
  sim_dat <- tryCatch(readRDS(data_filepath), error = function(e) {
    if (debug_mode) cat("DEBUG: Error loading RDS file:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(sim_dat) || !all(c("N", "T", "data") %in% names(sim_dat))) {
    return(paste("ERROR: Failed to load or invalid data structure in", fname_base))
  }
  N_val <- sim_dat$N
  T_val <- sim_dat$T
  df <- sim_dat$data

  # --- Prepare y_t list (vectorized format) ---
  y_t_list <- NULL
  data_prep_error_msg <- NULL
  tryCatch(
    {
      if (debug_mode) cat("DEBUG: Preparing y_t data structure...\n")
      df_sorted <- df %>% dplyr::arrange(i, t) # Use explicit dplyr::arrange
      expected_rows <- N_val * T_val
      if (nrow(df_sorted) != expected_rows) stop(sprintf("Data rows (%d) != N*T (%d).", nrow(df_sorted), expected_rows))
      if (length(unique(df_sorted$t)) != T_val || min(df_sorted$t) != 1 || max(df_sorted$t) != T_val) stop("Time points incorrect.")
      if (anyNA(df_sorted$y1) || anyNA(df_sorted$y2)) stop("NA values found in y1 or y2 columns.")
      if (!all(is.finite(df_sorted$y1)) || !all(is.finite(df_sorted$y2))) stop("Non-finite values (Inf/-Inf) found in y1 or y2 columns.")

      y_t_list <- vector("list", T_val)
      for (t_loop in 1:T_val) {
        y_matrix_t <- df_sorted %>%
          dplyr::filter(t == t_loop) %>%
          dplyr::select(y1, y2) %>%
          as.matrix()
        if (!all(dim(y_matrix_t) == c(N_val, 2))) stop(sprintf("Matrix slice dim error at t=%d. Expected %dx2, Got %s.", t_loop, N_val, paste(dim(y_matrix_t), collapse = "x")))
        if (any(!is.finite(y_matrix_t))) stop(sprintf("Non-finite values detected in matrix slice at t=%d.", t_loop)) # Extra check
        y_t_list[[t_loop]] <- y_matrix_t
      }
      if (length(y_t_list) != T_val) stop("Final list length incorrect.")
      if (debug_mode) cat("DEBUG: y_t data structure prepared successfully.\n")
    },
    error = function(e) {
      data_prep_error_msg <<- conditionMessage(e)
    }
  ) # Assign error message

  if (!is.null(data_prep_error_msg)) {
    return(paste("ERROR: Data prep failed for", fname_base, ":", data_prep_error_msg))
  }
  if (is.null(y_t_list)) {
    return(paste("ERROR: y_t_list is NULL after data prep for", fname_base))
  }

  # --- Stan data list ---
  stan_data <- list(N = N_val, T = T_val, y_t = y_t_list)
  if (is_multilevel) {
    # These priors are only expected by the ML models in their data block
    stan_data$prior_sigma_mu_scale <- prior_sigma_mu_scale
    stan_data$prior_sigma_phi_scale <- prior_sigma_phi_scale
  }

  if (debug_mode) {
    cat("DEBUG: Stan Data List Structure:\n")
    print(str(stan_data))
    cat("DEBUG: Summary of y_t[[2]] (first 5 rows):\n") # Look at an early time point matrix
    print(summary(stan_data$y_t[[2]][1:min(5, N_val), ]))
    cat("DEBUG: Any non-finite values in full y_t list? ", any(!is.finite(unlist(stan_data$y_t))), "\n")
    if (is_multilevel) {
      cat("DEBUG: ML prior scales passed: mu=", stan_data$prior_sigma_mu_scale, "phi=", stan_data$prior_sigma_phi_scale, "\n")
    } else {
      cat("DEBUG: SL model, no prior scales passed.\n")
    }
    cat("-----------------------------------------------------\n")
  }

  # --- Initialization Functions (Conditional based on is_multilevel) ---
  # Use stats:: namespace explicitly for clarity
  generate_inits_base <- function(N_val_init, p_mu_scale, p_phi_scale, is_ml) {
    base_list <- list(
      mu_global = stats::rnorm(2, 0, 0.1),
      rho = stats::runif(1, -0.5, 0.5)
    )
    if (is_ml) {
      base_list$phi11_base <- stats::runif(1, -0.5, 0.5)
      base_list$phi12_base <- stats::runif(1, -0.3, 0.3)
      base_list$phi21_base <- stats::runif(1, -0.3, 0.3)
      base_list$phi22_base <- stats::runif(1, -0.5, 0.5)
      base_list$sigma_mu <- abs(stats::rnorm(2, 0, p_mu_scale * 0.2)) + 0.01
      base_list$sigma_phi <- abs(stats::rnorm(4, 0, p_phi_scale * 0.2)) + 0.01
      # Ensure raw deviates are finite vectors of correct length N
      base_list$mu1_raw <- stats::rnorm(N_val_init, 0, 1)
      base_list$mu2_raw <- stats::rnorm(N_val_init, 0, 1)
      base_list$dev_phi11_raw <- stats::rnorm(N_val_init, 0, 1)
      base_list$dev_phi12_raw <- stats::rnorm(N_val_init, 0, 1)
      base_list$dev_phi21_raw <- stats::rnorm(N_val_init, 0, 1)
      base_list$dev_phi22_raw <- stats::rnorm(N_val_init, 0, 1)
    } else {
      base_list$phi11 <- stats::runif(1, -0.5, 0.5)
      base_list$phi12 <- stats::runif(1, -0.3, 0.3)
      base_list$phi21 <- stats::runif(1, -0.3, 0.3)
      base_list$phi22 <- stats::runif(1, -0.5, 0.5)
    }
    # Check for non-finite values in the generated list
    if (any(!sapply(base_list, function(x) all(is.finite(x))))) {
      stop("Non-finite values generated in base_list during initialization.")
    }
    return(base_list)
  }
  init_fun_norm <- function(chain_id = 1, common_inits_base, data_df) {
    sd_y1 <- max(0.1, stats::sd(data_df$y1, na.rm = TRUE))
    sd_y2 <- max(0.1, stats::sd(data_df$y2, na.rm = TRUE))
    init_sigma1 <- sd_y1 * stats::runif(1, 0.5, 1.5)
    init_sigma2 <- sd_y2 * stats::runif(1, 0.5, 1.5)
    # Ensure sigma is positive
    if (init_sigma1 <= 0 || init_sigma2 <= 0) stop("Generated sigma <= 0 during initialization.")
    specific_inits <- list(sigma = c(init_sigma1, init_sigma2))
    full_init <- c(common_inits_base, specific_inits)
    # Final check
    if (any(!sapply(full_init, function(x) all(is.finite(x))))) {
      stop("Non-finite values generated in init_fun_norm.")
    }
    return(full_init)
  }
  init_fun_skewnorm <- function(chain_id = 1, common_inits_base, data_df) {
    sd_y1 <- max(0.1, stats::sd(data_df$y1, na.rm = TRUE))
    sd_y2 <- max(0.1, stats::sd(data_df$y2, na.rm = TRUE))
    # Ensure omega > 0
    init_omega1 <- max(0.01, sd_y1 * stats::runif(1, 0.5, 1.5))
    init_omega2 <- max(0.01, sd_y2 * stats::runif(1, 0.5, 1.5))
    # *** Tighter Inits for alpha and xi ***
    init_alpha1 <- stats::rnorm(1, 0, 1.0) # Reduced SD from 2.0 to 1.0
    init_alpha2 <- stats::rnorm(1, 0, 1.0) # Reduced SD from 2.0 to 1.0
    init_xi1 <- stats::rnorm(1, 0, 0.1) # Reduced SD from 0.2 to 0.1
    init_xi2 <- stats::rnorm(1, 0, 0.1) # Reduced SD from 0.2 to 0.1

    specific_inits <- list(omega = c(init_omega1, init_omega2), alpha = c(init_alpha1, init_alpha2), xi = c(init_xi1, init_xi2))
    full_init <- c(common_inits_base, specific_inits)
    # Final check for finite values
    if (any(!sapply(full_init, function(x) all(is.finite(x))))) {
      stop("Non-finite values generated in init_fun_skewnorm.")
    }
    return(full_init)
  }

  # --- Generate Inits ---
  init_list_norm_vec <- NULL
  init_list_skewnorm <- NULL
  init_error_msg <- NULL
  tryCatch(
    {
      if (debug_mode) cat("DEBUG: Generating initial values...\n")
      # Pass the is_multilevel flag to the base init function
      base_inits_list <- lapply(1:stan_chains, function(id) generate_inits_base(N_val, prior_sigma_mu_scale, prior_sigma_phi_scale, is_multilevel))
      init_list_norm_vec <- lapply(1:stan_chains, function(id) init_fun_norm(id, base_inits_list[[id]], df))
      init_list_skewnorm <- lapply(1:stan_chains, function(id) init_fun_skewnorm(id, base_inits_list[[id]], df))
      if (debug_mode) {
        cat("DEBUG: Initial values generated successfully.\n")
        cat("DEBUG: Init list for SkewNormal (Chain 1):\n")
        print(init_list_skewnorm[[1]])
        cat("DEBUG: Init list for Normal (Chain 1):\n")
        print(init_list_norm_vec[[1]])
        cat("-----------------------------------------------------\n")
      }
    },
    error = function(e) {
      init_error_msg <<- conditionMessage(e)
    }
  ) # Assign NULL on error

  if (!is.null(init_error_msg)) {
    return(paste("ERROR: Init generation failed for", fname_base, ":", init_error_msg))
  }
  if (is.null(init_list_norm_vec) || is.null(init_list_skewnorm)) {
    return(paste("ERROR: Init lists are NULL after generation for", fname_base))
  }

  # --- Fit Models ---
  stan_control <- list(adapt_delta = stan_adapt_delta, max_treedepth = stan_max_treedepth)
  common_seed_base <- 2024 * cond_num + rep_num # Seed based on file IDs
  fit_results <- list(norm_status = "Not run", skew_status = "Not run") # Store results for final message

  # --- Fit Skew-Normal Model ---
  if (!file.exists(out_file_sn)) {
    if (debug_mode) cat(sprintf("DEBUG: Attempting SkewNormal fit (Seed: %d)...\n", common_seed_base + 1000))
    fit_sn <- NULL
    skew_error_msg <- NULL
    start_time <- Sys.time()
    tryCatch(
      {
        # Use suppressMessages to reduce console clutter from Stan iterations unless debugging
        # Use refresh=0 if not debugging and running parallel to avoid interleaved messages
        current_refresh <- if (debug_mode) 500 else 0
        fit_sn <- rstan::sampling(
          object = compiled_skew_model, data = stan_data,
          chains = stan_chains, iter = stan_iter, warmup = stan_warmup,
          seed = common_seed_base + 1000, init = init_list_skewnorm,
          refresh = current_refresh, control = stan_control
        )
        if (debug_mode) cat("DEBUG: SkewNormal sampling call finished.\n")
        # Save model fit object immediately
        saveRDS(fit_sn, file = out_file_sn)
        if (debug_mode) cat("DEBUG: SkewNormal fit saved to:", out_file_sn, "\n")
      },
      error = function(e) {
        # Capture the specific error from Stan
        skew_error_msg <<- paste("Stan Sampling Error:", conditionMessage(e))
        if (debug_mode) {
          cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
          cat("DEBUG: ERROR during SkewNormal sampling:\n")
          print(e) # Print the full error object
          cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        }
      }
    ) # End tryCatch
    end_time <- Sys.time()
    elapsed_time <- difftime(end_time, start_time, units = "mins")

    if (!is.null(skew_error_msg)) {
      fit_results$skew_status <- skew_error_msg # Use the captured message
    } else {
      # Check samples only if no error occurred during sampling call
      # Need to handle cases where sampling "finishes" but is empty (e.g., init failure within C++)
      samples_ok <- FALSE
      if (!is.null(fit_sn)) {
        # Use summary() to check for major issues post-sampling
        summary_sn <- tryCatch(summary(fit_sn)$summary, error = function(e) NULL)
        if (!is.null(summary_sn) && nrow(summary_sn) > 0 && all(is.finite(summary_sn[, "mean"]))) {
          samples_ok <- TRUE
        }
      }

      fit_results$skew_status <- if (samples_ok) {
        sprintf("SkewNormal OK (%.1f min)", elapsed_time)
      } else {
        "SkewNormal sampling FAILED (No valid samples/summary)"
      }
    }
  } else {
    fit_results$skew_status <- "SkewNormal skipped (exists)"
    if (debug_mode) cat("DEBUG: SkewNormal fit skipped (file exists).\n")
  }
  if (debug_mode) cat("-----------------------------------------------------\n")

  # --- Fit Normal Model ---
  if (!file.exists(out_file_norm)) {
    if (debug_mode) cat(sprintf("DEBUG: Attempting Normal fit (Seed: %d)...\n", common_seed_base + 2000))
    fit_norm <- NULL
    norm_error_msg <- NULL
    start_time <- Sys.time()
    tryCatch(
      {
        current_refresh <- if (debug_mode) 500 else 0
        fit_norm <- rstan::sampling(
          object = compiled_norm_model, data = stan_data,
          chains = stan_chains, iter = stan_iter, warmup = stan_warmup,
          seed = common_seed_base + 2000, init = init_list_norm_vec,
          refresh = current_refresh, control = stan_control
        )
        if (debug_mode) cat("DEBUG: Normal sampling call finished.\n")
        saveRDS(fit_norm, file = out_file_norm)
        if (debug_mode) cat("DEBUG: Normal fit saved to:", out_file_norm, "\n")
      },
      error = function(e) {
        norm_error_msg <<- paste("Stan Sampling Error:", conditionMessage(e))
        if (debug_mode) {
          cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
          cat("DEBUG: ERROR during Normal sampling:\n")
          print(e) # Print the full error object
          cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        }
      }
    ) # End tryCatch
    end_time <- Sys.time()
    elapsed_time <- difftime(end_time, start_time, units = "mins")

    if (!is.null(norm_error_msg)) {
      fit_results$norm_status <- norm_error_msg
    } else {
      samples_ok_norm <- FALSE
      if (!is.null(fit_norm)) {
        summary_norm <- tryCatch(summary(fit_norm)$summary, error = function(e) NULL)
        if (!is.null(summary_norm) && nrow(summary_norm) > 0 && all(is.finite(summary_norm[, "mean"]))) {
          samples_ok_norm <- TRUE
        }
      }
      fit_results$norm_status <- if (samples_ok_norm) {
        sprintf("Normal OK (%.1f min)", elapsed_time)
      } else {
        "Normal sampling FAILED (No valid samples/summary)"
      }
    }
  } else {
    fit_results$norm_status <- "Normal skipped (exists)"
    if (debug_mode) cat("DEBUG: Normal fit skipped (file exists).\n")
  }

  final_message <- paste("Processed:", fname_base, "-", fit_results$norm_status, ";", fit_results$skew_status)
  if (debug_mode) {
    cat("-----------------------------------------------------\n")
    cat("DEBUG: Final Status:", final_message, "\n")
    cat("=====================================================\n\n")
  }
  # Return the combined status message
  return(final_message)
} # End of fit_models_for_single_file_worker


# Ensure necessary libraries are loaded before calling this function
# Requires: rstan, dplyr, tidyr, foreach, doParallel, parallel, stats

# Assumes the helper function fit_models_for_single_file_worker is defined
# previously in the script or environment, matching the version from
# https://chat.google.com/share/UnGkMqkH6c3H (or similar) which includes
# the debug_mode argument and detailed print statements.

# --- Main Fitting Function (Called by run_pipeline.R) ---

#' Fit Normal and Skew-Normal models to simulated data (Parallelized Internally)
#' Handles selection between multilevel and single-level models based on condition.
#' Includes debug mode for sequential file processing WITH parallel chains,
#' or normal mode for parallel file processing with sequential chains.
#'
#' @param data_dir Directory containing the simulation RDS files (e.g., data_sn).
#' @param fits_dir Directory where the Stan fit RDS files will be saved (e.g., fits_sn).
#' @param stan_models_dir Directory where the .stan model files are located.
#' @param prior_sigma_mu_scale Scale for the prior on sigma_mu (for ML models).
#' @param prior_sigma_phi_scale Scale for the prior on sigma_phi (for ML models).
#' @param stan_iter Total iterations per chain.
#' @param stan_warmup Warmup iterations per chain.
#' @param stan_chains Number of chains per model fit. Also used to set mc.cores in debug mode.
#' @param stan_adapt_delta Target acceptance probability for NUTS.
#' @param stan_max_treedepth Maximum tree depth for NUTS.
#' @param num_cores Number of cores for file-level parallelism (used only when debug_mode=FALSE).
#' @param debug_mode Logical. If TRUE, run sequentially across files BUT use parallel chains within Stan fits.
#' @param debug_file_limit Integer. Maximum number of files to process if debug_mode=TRUE. NULL to process all.
#' @return Invisibly returns NULL. Saves fit files to disk and prints progress.
fit_multilevel_models_sn_vec <- function(
    data_dir,
    fits_dir,
    stan_models_dir,
    prior_sigma_mu_scale = 1.0,
    prior_sigma_phi_scale = 0.5,
    stan_iter = 3000,
    stan_warmup = 1500,
    stan_chains = 4, # Number of chains Stan will run
    stan_adapt_delta = 0.80,
    stan_max_treedepth = 10,
    num_cores = 8, # Cores for file-level parallelism (used only when debug_mode=FALSE)
    debug_mode = FALSE,
    debug_file_limit = 5) {
  # --- Helper functions ---
  # Null default operator
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # --- Set Stan Options ---
  rstan::rstan_options(auto_write = TRUE)

  # *** Manage mc.cores based on debug_mode ***
  original_mc_cores <- getOption("mc.cores", 1L) # Store original global setting
  cores_to_set <- 1 # Default if not debugging

  if (debug_mode) {
    # DEBUG MODE: Set mc.cores to run chains in parallel.
    # Use the number of chains requested for the model fit,
    # but don't exceed the number of physical cores available.
    # Use logical=FALSE for physical cores if available.
    # Ensure parallel package is available for detectCores
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' needed to detect cores for mc.cores setting. Defaulting to stan_chains.")
      available_cores <- stan_chains
    } else {
      available_cores <- parallel::detectCores(logical = FALSE) %||% 1L
    }
    cores_for_chains <- min(stan_chains, available_cores)
    # Ensure at least 1 core
    cores_to_set <- max(1L, cores_for_chains)
    cat(sprintf("--- DEBUG MODE: Setting options(mc.cores = %d) to run chains in parallel for each fit. ---\n", cores_to_set))
  } else {
    # NORMAL (non-debug) MODE: Set mc.cores to 1.
    # Parallelism is handled by foreach across files.
    cores_to_set <- 1
    cat(sprintf("--- NORMAL MODE: Setting options(mc.cores = %d). Parallelism via foreach across files enabled. ---\n", cores_to_set))
  }
  options(mc.cores = cores_to_set)
  # Ensure original setting is restored when function exits
  on.exit(options(mc.cores = original_mc_cores), add = TRUE) # Add=TRUE ensures it doesn't overwrite other on.exit calls

  # --- Validate Directories ---
  if (!dir.exists(data_dir)) stop("Data directory not found: ", data_dir)
  if (!dir.exists(stan_models_dir)) stop("Stan models directory not found: ", stan_models_dir)
  if (!dir.exists(fits_dir)) {
    warning("Fits directory not found, creating: ", fits_dir)
    tryCatch(dir.create(fits_dir, recursive = TRUE),
      error = function(e) stop("Failed to create fits directory: ", conditionMessage(e))
    )
  }

  # --- Stan Model Filenames (ALL FOUR) ---
  stan_ml_norm_file <- file.path(stan_models_dir, "model_normal_multilevel_vec.stan")
  stan_ml_skew_file <- file.path(stan_models_dir, "model_skewnormal_multilevel_vec.stan")
  stan_sl_norm_file <- file.path(stan_models_dir, "model_normal_singlelevel_vec.stan")
  stan_sl_skew_file <- file.path(stan_models_dir, "model_skewnormal_singlelevel_vec.stan")
  # Check existence
  if (!file.exists(stan_ml_norm_file)) stop("Stan file not found: ", stan_ml_norm_file)
  if (!file.exists(stan_ml_skew_file)) stop("Stan file not found: ", stan_ml_skew_file)
  if (!file.exists(stan_sl_norm_file)) stop("Stan file not found: ", stan_sl_norm_file)
  if (!file.exists(stan_sl_skew_file)) stop("Stan file not found: ", stan_sl_skew_file)


  # --- Compile ALL Models ONCE ---
  # Consider adding caching here if compilation is slow and models don't change often
  cat("Compiling Stan models (ML & SL, if needed)...\n")
  sm_norm_ml <- tryCatch(rstan::stan_model(stan_ml_norm_file, model_name = "ML_Norm_Vec"), error = function(e) stop("Compile ML Norm failed:", e$message))
  sm_skew_ml <- tryCatch(rstan::stan_model(stan_ml_skew_file, model_name = "ML_Skew_Vec"), error = function(e) stop("Compile ML Skew failed:", e$message))
  sm_norm_sl <- tryCatch(rstan::stan_model(stan_sl_norm_file, model_name = "SL_Norm_Vec"), error = function(e) stop("Compile SL Norm failed:", e$message))
  sm_skew_sl <- tryCatch(rstan::stan_model(stan_sl_skew_file, model_name = "SL_Skew_Vec"), error = function(e) stop("Compile SL Skew failed:", e$message))
  cat("Stan models compiled.\n")

  # --- Load Simulation Conditions (for lookup) ---
  # Construct path assuming data_dir is like ".../data" and conditions are in ".../data/sim_conditions.rds"
  # Adjust if your structure is different
  conditions_file <- file.path(dirname(data_dir), "data", "sim_conditions.rds") # Assumes data_dir is one level down
  if (!file.exists(conditions_file)) {
    # Fallback if data_dir *is* the base data directory
    conditions_file <- file.path(data_dir, "sim_conditions.rds")
    if (!file.exists(conditions_file)) {
      stop(paste("sim_conditions.rds not found in", data_dir, "or parent's data directory."))
    }
  }
  sim_conditions_lookup <- tryCatch(readRDS(conditions_file),
    error = function(e) stop("Failed to read sim_conditions.rds: ", conditionMessage(e))
  )

  # Basic validation of lookup table
  if (!all(c("condition_id", "has_random_effects") %in% names(sim_conditions_lookup))) {
    stop("sim_conditions.rds must contain 'condition_id' and 'has_random_effects' columns.")
  }
  sim_conditions_lookup <- sim_conditions_lookup %>%
    dplyr::select(condition_id, has_random_effects) %>%
    dplyr::distinct()

  # --- Get List of Data Files ---
  all_data_files <- list.files(data_dir, pattern = "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (length(all_data_files) == 0) {
    cat("No sim_data files found in: ", data_dir, "\n")
    return(invisible(NULL))
  }

  # --- Apply Debug Limit ---
  if (debug_mode) {
    limit_msg <- "(Processing all files sequentially - debug_file_limit is NULL or <= 0)"
    files_to_keep <- length(all_data_files)
    if (!is.null(debug_file_limit) && debug_file_limit > 0) {
      files_to_keep <- min(debug_file_limit, length(all_data_files))
      limit_msg <- sprintf("(Processing only the first %d file(s) sequentially - foreach %%do%%)", files_to_keep)
    }
    data_files_to_process <- all_data_files[1:files_to_keep]
    cat(sprintf("!!! DEBUG MODE ACTIVE: %s !!!\n", limit_msg))
  } else {
    data_files_to_process <- all_data_files
    cat(sprintf("Found %d simulation data files to process.\n", length(data_files_to_process)))
  }
  if (length(data_files_to_process) == 0) {
    cat("No data files selected for processing.\n")
    return(invisible(NULL))
  }

  # --- Setup Parallel Backend (or run sequentially in debug mode) ---
  # Cluster variable needs to be defined for cleanup even if not used in debug mode
  cl <- NULL
  run_engine <- NULL # Define execution engine

  if (!debug_mode) {
    # NORMAL MODE: Setup parallel backend for foreach %dopar%
    n_cores <- max(1, num_cores) # Use the num_cores argument for file-level parallelism
    cat(sprintf("Setting up parallel cluster with %d cores for file processing...\n", n_cores))
    # Check if necessary packages are installed
    if (!requireNamespace("doParallel", quietly = TRUE) || !requireNamespace("foreach", quietly = TRUE)) {
      stop("Packages 'doParallel' and 'foreach' are required for parallel execution (when debug_mode=FALSE).")
    }

    # Check if a cluster is already registered, stop if needed
    if (foreach::getDoParRegistered()) {
      cat("Warning: Existing parallel cluster detected. Stopping it before creating a new one.\n")
      # Try to stop cluster associated with 'cl' if it exists in the environment
      if (exists("cl_prev") && inherits(cl_prev, "cluster")) try(parallel::stopCluster(cl_prev), silent = TRUE)
      try(doParallel::stopImplicitCluster(), silent = TRUE) # Stop any implicit cluster
    }
    # Create and register new cluster
    cl <- tryCatch(parallel::makeCluster(n_cores),
      error = function(e) stop("Failed to create parallel cluster: ", conditionMessage(e))
    )
    doParallel::registerDoParallel(cl)
    cat("Cluster registered.\n")

    # Ensure cluster is stopped when the function exits (normally or on error)
    on.exit(
      {
        if (!is.null(cl)) {
          cat("Stopping parallel cluster...\n")
          try(parallel::stopCluster(cl), silent = TRUE)
          try(doParallel::registerDoSEQ(), silent = TRUE) # Deregister parallel backend
          cat("Parallel cluster stopped.\n")
        }
        # mc.cores is restored by the separate on.exit call added earlier
      },
      add = TRUE
    ) # Add this cleanup action

    run_engine <- foreach::`%dopar%` # Use parallel engine
  } else {
    # DEBUG MODE: Use sequential foreach %do%
    cat("Running sequentially (foreach %do%) in debug mode.\n")
    run_engine <- foreach::`%do%` # Use sequential engine
    n_cores <- 1 # No file-level parallelism
  }


  # --- Run Loop (Parallel or Sequential based on debug_mode) ---
  cat(sprintf("Starting fitting (%s)...\n", if (debug_mode) "Sequential Files (Parallel Chains)" else "Parallel Files (Sequential Chains)"))
  start_time_loop <- Sys.time()

  # Explicitly list ALL variables/functions needed by the worker that are defined *outside* it
  # This is crucial for parallel execution (%dopar%)
  # Make sure fit_models_for_single_file_worker is available in the execution environment
  # If it's defined globally before calling this main function, it should be found.
  # If defined locally, it needs to be exported. Assuming global/sourced definition here.
  if (!exists("fit_models_for_single_file_worker", mode = "function")) {
    stop("Helper function 'fit_models_for_single_file_worker' not found. Ensure it is defined or sourced before calling this function.")
  }
  vars_to_export <- c( # Keep this minimal if worker function is global
    "sm_norm_ml", "sm_skew_ml", "sm_norm_sl", "sm_skew_sl",
    "sim_conditions_lookup", "fits_dir",
    "prior_sigma_mu_scale", "prior_sigma_phi_scale",
    "stan_iter", "stan_warmup", "stan_chains",
    "stan_adapt_delta", "stan_max_treedepth",
    "debug_mode"
  )
  # If fit_models_for_single_file_worker IS NOT defined globally, add its name here:
  vars_to_export <- c("fit_models_for_single_file_worker", vars_to_export)


  # Prepare foreach object
  # Ensure necessary packages for the WORKER function are listed
  foreach_obj <- foreach::foreach(
    f = data_files_to_process,
    .packages = c("rstan", "dplyr", "tidyr", "stats"), # Worker needs these
    .export = vars_to_export,
    .errorhandling = "pass", # Return error object if worker fails
    .verbose = debug_mode # Show some foreach progress if debugging
  )

  # Execute using the chosen engine (%do% or %dopar%)
  results_list <- run_engine(foreach_obj, {
    # --- Code executed by each worker (or sequentially) ---

    # Need to load the worker function if it wasn't defined globally
    # Assuming it's available via export or global environment here

    # Parse condition ID from filename 'f'
    fname_base_par <- basename(f)
    cond_match_par <- regexpr("(?<=_cond)(\\d+)", fname_base_par, perl = TRUE)
    if (cond_match_par == -1 || attr(cond_match_par, "match.length") == -1) {
      # Stop this specific worker/iteration if parsing fails
      # Return an error object instead of stopping the whole cluster
      return(simpleError(paste0("FATAL: Cannot parse condition ID in worker for file: ", fname_base_par)))
    }
    cond_num_par <- tryCatch(as.integer(regmatches(fname_base_par, cond_match_par)), error = function(e) NA)
    if (is.na(cond_num_par)) {
      return(simpleError(paste0("FATAL: Cannot convert condition ID to integer in worker for file: ", fname_base_par)))
    }

    # Look up the flag from the exported lookup table
    # Use tryCatch for robustness in case lookup fails
    cond_row <- tryCatch(sim_conditions_lookup[sim_conditions_lookup$condition_id == cond_num_par, ],
      error = function(e) NULL
    )
    if (is.null(cond_row) || nrow(cond_row) != 1) {
      return(simpleError(paste0("FATAL: Condition ID ", cond_num_par, " not found or duplicated in lookup table.")))
    }
    is_ml_flag <- cond_row$has_random_effects[1]

    # Select the correct pair of compiled models to pass
    compiled_norm_to_use <- if (is_ml_flag) sm_norm_ml else sm_norm_sl
    compiled_skew_to_use <- if (is_ml_flag) sm_skew_ml else sm_skew_sl

    # Call the worker function
    # rstan::sampling inside the worker will respect the mc.cores option set earlier by the main function
    # The worker function must be available in this scope
    fit_models_for_single_file_worker(
      data_filepath = f,
      compiled_norm_model = compiled_norm_to_use,
      compiled_skew_model = compiled_skew_to_use,
      fits_dir = fits_dir,
      is_multilevel = is_ml_flag,
      prior_sigma_mu_scale = prior_sigma_mu_scale,
      prior_sigma_phi_scale = prior_sigma_phi_scale,
      stan_iter = stan_iter,
      stan_warmup = stan_warmup,
      stan_chains = stan_chains,
      stan_adapt_delta = stan_adapt_delta,
      stan_max_treedepth = stan_max_treedepth,
      debug_mode = debug_mode # Pass the flag through
    )
  }) # End foreach loop execution

  end_time_loop <- Sys.time()
  cat(sprintf(
    "\nLoop finished. Total time: %.2f minutes.\n",
    difftime(end_time_loop, start_time_loop, units = "mins")
  ))

  # --- Process Results (Summarize success/failure) ---
  cat("Processing fitting results...\n")
  errors <- list()
  success_count <- 0
  skipped_count <- 0
  status_messages <- character(length(results_list)) # Store status messages

  for (i in seq_along(results_list)) {
    res <- results_list[[i]]
    current_file <- data_files_to_process[i] # Get corresponding filename

    if (inherits(res, "error")) {
      # This catches errors from the foreach worker itself (e.g., parsing file, lookup, or FATAL errors returned)
      err_msg <- paste("FATAL Worker Error:", conditionMessage(res))
      errors[[length(errors) + 1]] <- list(file = basename(current_file), error = err_msg)
      status_messages[i] <- err_msg
    } else if (is.character(res)) {
      # This is the status string returned by fit_models_for_single_file_worker
      status_messages[i] <- res
      # Parse the status message for summary counts and errors reported *by the worker*
      if (grepl("^Skipped:", res)) skipped_count <- skipped_count + 1
      # Count file as successfully processed if no fatal error and at least one model OK or skipped
      if (!grepl("FATAL", res) && (grepl(" OK", res) || grepl("skipped", res))) {
        success_count <- success_count + 1
      }
      # Capture non-fatal errors reported within the worker's return string
      if (grepl("ERROR:", res) || grepl("FAILED", res) || grepl("Error:", res)) { # More robust check
        # Avoid double-counting fatal errors already caught above
        is_fatal <- any(sapply(errors, function(e) e$file == basename(current_file) && grepl("FATAL", e$error)))
        if (!is_fatal) {
          errors[[length(errors) + 1]] <- list(file = basename(current_file), error = res)
        }
      }
    } else {
      # Should not happen if worker returns string or foreach passes error
      err_msg <- paste("Unknown result type from worker for file:", basename(current_file))
      errors[[length(errors) + 1]] <- list(file = basename(current_file), error = err_msg)
      status_messages[i] <- err_msg
    }
  }

  # Print all status messages (can be long)
  cat("--- Individual File Status ---\n")
  if (length(status_messages) > 0) {
    for (msg in status_messages) {
      cat("   ", msg, "\n")
    }
  } else {
    cat("   No results messages generated.\n")
  }
  cat("-----------------------------\n")

  total_processed_or_skipped <- length(data_files_to_process) # Count based on files attempted
  cat(sprintf(
    "Fitting Summary: %d files attempted, %d reported errors.\n",
    total_processed_or_skipped, length(errors)
  ))

  if (length(errors) > 0) {
    cat("Errors encountered during fitting:\n")
    # Deduplicate errors based on file, keeping the most severe ("FATAL" > "ERROR" > "FAILED")
    error_summary <- list()
    for (err in errors) {
      fname <- err$file
      current_msg <- err$error
      severity <- ifelse(grepl("FATAL", current_msg), 3, ifelse(grepl("ERROR", current_msg) | grepl("Error:", current_msg), 2, ifelse(grepl("FAILED", current_msg), 1, 0)))
      if (!fname %in% names(error_summary) || severity > error_summary[[fname]]$severity) {
        error_summary[[fname]] <- list(msg = current_msg, severity = severity)
      }
    }
    # Print summarized errors
    max_errors_to_print <- 50
    printed_count <- 0
    # Create data frame for sorting
    if (length(error_summary) > 0) {
      error_df <- data.frame(
        file = names(error_summary),
        msg = sapply(error_summary, `[[`, "msg"),
        severity = sapply(error_summary, `[[`, "severity"),
        stringsAsFactors = FALSE
      )
      error_df <- error_df[order(-error_df$severity, error_df$file), ] # Sort by severity, then name

      for (i in 1:nrow(error_df)) {
        if (printed_count >= max_errors_to_print) break
        cat(sprintf("   - File: %s, Detail: %s\n", error_df$file[i], error_df$msg[i]))
        printed_count <- printed_count + 1
      }
      if (nrow(error_df) > max_errors_to_print) cat(sprintf("   (... %d more errors not printed ...)\n", nrow(error_df) - max_errors_to_print))
    } else {
      cat("   (No unique errors found in summary - check status messages)\n")
    }
  }

  # The on.exit handlers will restore the original mc.cores setting and stop the cluster if needed.
  cat("\nFitting function completed.\n")
  invisible(NULL) # Return nothing
} # End of main function fit_multilevel_models_sn_vec
