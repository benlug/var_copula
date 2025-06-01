###########################################################################
# process_fits.R
#
# Processes Stan fits from VAR(1) Copula simulation (handles separate alphas).
# Reads fits from ../fits, data/conditions from ../data.
# Saves aggregated summaries to the current 'results' directory.
###########################################################################

# --- Libraries ---
library(rstan)
library(tidyverse)
library(posterior) # Potentially useful for more advanced summaries later
library(this.path)
library(stringr)

# --- Configuration ---
cat("--- Starting Fit Processing ---\n")
# Determine script's directory to set relative paths correctly
# Use getwd() if running interactively from the 'results' directory
# Use this.dir() if sourcing the script via run_pipeline.R
RESULTS_DIR <- tryCatch(this.dir(), error = function(e) getwd())
setwd(RESULTS_DIR) # Ensure working directory is results for relative paths

FITS_DIR <- file.path(RESULTS_DIR, "../fits")
DATA_DIR <- file.path(RESULTS_DIR, "../data")
if (!dir.exists(FITS_DIR)) stop("Fits directory not found: ", normalizePath(FITS_DIR))
if (!dir.exists(DATA_DIR)) stop("Data directory not found: ", normalizePath(DATA_DIR))

# --- Load Simulation Conditions ---
sim_conds_file <- file.path(DATA_DIR, "sim_conditions.rds") # Simplified name
if (!file.exists(sim_conds_file)) stop("Sim conditions file not found: ", sim_conds_file)
sim_conditions_df <- readRDS(sim_conds_file) %>% mutate(condition_id = as.integer(condition_id))
cat("Loaded", nrow(sim_conditions_df), "simulation conditions.\n")

# --- Helper Functions ---
`%||%` <- function(a, b) if (!is.null(a)) a else b

# modified helper that also returns reason when a fit can't be read
safe_read_stanfit <- function(filename) {
  if (!file.exists(filename)) {
    return(list(fit = NULL, reason = "file missing"))
  }
  fit_obj <- tryCatch(readRDS(filename), error = function(e) e)
  if (inherits(fit_obj, "error")) {
    return(list(fit = NULL, reason = paste("readRDS error:", fit_obj$message)))
  }
  if (!inherits(fit_obj, "stanfit")) {
    cls <- paste(class(fit_obj), collapse = ",")
    return(list(fit = NULL, reason = paste("not stanfit (", cls, ")")))
  }
  samples_ok <- FALSE
  tryCatch({
    sims <- slot(fit_obj, "sim")$samples
    if (is.list(sims) && length(sims) > 0 && length(sims[[1]]) > 0) {
      samples_ok <- TRUE
    }
  }, error = function(e) {})
  if (!samples_ok) {
    return(list(fit = NULL, reason = "stanfit has no samples"))
  }
  list(fit = fit_obj, reason = NULL)
}

calc_e_fmi <- function(energy_vec) {
  # compute the energy fraction of missing information
  n <- length(energy_vec)
  if (n < 3) {
    return(NA_real_)
  }
  energy_vec <- energy_vec[is.finite(energy_vec)]
  n <- length(energy_vec)
  if (n < 3) {
    return(NA_real_)
  }
  Emean <- mean(energy_vec)
  varE <- stats::var(energy_vec)
  if (n < 2 || is.na(varE) || varE <= 0) {
    return(NA_real_)
  } # Added checks for varE
  lag1_sum <- sum((energy_vec[-1] - Emean) * (energy_vec[-n] - Emean))
  if (!is.finite(lag1_sum)) {
    return(NA_real_)
  }
  lag1 <- lag1_sum / (n - 1)
  denom <- varE + 2 * lag1
  if (is.na(denom) || denom <= 1e-9) {
    return(NA_real_)
  }
  efmi_val <- varE / denom
  pmax(0, pmin(1.5, efmi_val))
}

get_true_value_var1 <- function(param_name, true_params_list, fitted_model_type) {
  # return the true parameter value for comparison
  true_val <- NA_real_
  true_margin1_dist <- true_params_list$margin1_dist %||% NA_character_
  true_margin2_dist <- true_params_list$margin2_dist %||% NA_character_
  true_margin1_params <- true_params_list$margin1_params
  true_margin2_params <- true_params_list$margin2_params

  # Added check if params lists exist
  if (is.null(true_margin1_params) || is.null(true_margin2_params)) {
    return(NA_real_)
  }

  if (param_name == "phi11") {
    true_val <- true_params_list$phi[1, 1]
  } else if (param_name == "phi12") {
    true_val <- true_params_list$phi[1, 2]
  } else if (param_name == "phi21") {
    true_val <- true_params_list$phi[2, 1]
  } else if (param_name == "phi22") {
    true_val <- true_params_list$phi[2, 2]
  } else if (param_name == "mu[1]") {
    true_val <- true_params_list$mu[1]
  } else if (param_name == "mu[2]") {
    true_val <- true_params_list$mu[2]
  } else if (grepl("^(sigma|xi|omega|alpha)\\[1\\]$", param_name)) {
    param_base <- sub("\\[1\\]", "", param_name)
    if (param_base == "sigma") {
      true_val <- ifelse(true_margin1_dist == "normal", true_margin1_params$sd, NA_real_)
    } else if (true_margin1_dist == "skewnormal") {
      true_val <- true_margin1_params[[param_base]] %||% NA_real_
    }
  } else if (grepl("^(sigma|xi|omega|alpha)\\[2\\]$", param_name)) {
    param_base <- sub("\\[2\\]", "", param_name)
    if (param_base == "sigma") {
      true_val <- ifelse(true_margin2_dist == "normal", true_margin2_params$sd, NA_real_)
    } else if (true_margin2_dist == "skewnormal") {
      true_val <- true_margin2_params[[param_base]] %||% NA_real_
    }
  } else if (param_name == "rho") {
    true_val <- ifelse(true_params_list$copula_type == "gaussian", true_params_list$copula_param_value, NA_real_)
  } else if (param_name == "theta") {
    true_val <- ifelse(true_params_list$copula_type == "clayton", true_params_list$copula_param_value, NA_real_)
  } else if (param_name == "tau") {
    true_val <- true_params_list$copula_tau
  }
  return(true_val)
}

get_param_category_var1 <- function(param_name) {
  # group parameters for summarizing
  if (grepl("^phi\\d{2}$", param_name)) {
    "VAR Coeffs"
  } else if (grepl("^mu\\[\\d\\]$", param_name)) {
    "Intercepts"
  } else if (grepl("^(sigma|xi|omega|alpha)\\[\\d\\]$", param_name)) {
    "Residual Params"
  } else if (grepl("^(rho|theta|tau)$", param_name)) {
    "Copula Params"
  } else {
    "Other"
  }
}

# --- Find Fit Files ---
fit_files <- list.files(FITS_DIR, pattern = "^fit_[A-Z]{2}_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
if (length(fit_files) == 0) stop("No fit files found in: ", normalizePath(FITS_DIR))
cat("Found", length(fit_files), "potential fit files.\n")

# --- Process Fits ---
param_results_list <- list()
sampler_results_list <- list()
files_processed <- 0
files_failed <- 0
failure_logs <- list()

for (fit_file in fit_files) {
  # cat(".") # Removed progress dot
  base_name <- basename(fit_file)
  match_info <- stringr::str_match(base_name, "^fit_([A-Z]{2})_cond(\\d+)_rep(\\d+)\\.rds$")
  if (is.na(match_info[1, 1])) {
    files_failed <- files_failed + 1
    next
  }
  fitted_model_code <- match_info[1, 2]
  cond_id <- as.integer(match_info[1, 3])
  rep_i <- as.integer(match_info[1, 4])

  fit_res <- safe_read_stanfit(fit_file)
  if (is.null(fit_res$fit)) {
    files_failed <- files_failed + 1
    failure_logs[[length(failure_logs) + 1]] <- tibble(
      file = base_name,
      condition_id = cond_id,
      rep_i = rep_i,
      reason = fit_res$reason
    )
    next
  }
  fit <- fit_res$fit

  sim_data_file <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cond_id, rep_i)) # Simplified name
  if (!file.exists(sim_data_file)) {
    files_failed <- files_failed + 1
    failure_logs[[length(failure_logs) + 1]] <- tibble(
      file = base_name,
      condition_id = cond_id,
      rep_i = rep_i,
      reason = "sim data missing"
    )
    next
  }
  sim_data <- readRDS(sim_data_file)
  true_params <- sim_data$true_params
  dgp_info <- sim_data$dgp_info
  # Added check for true_params structure
  if (is.null(true_params) || is.null(dgp_info) || is.null(true_params$margin1_params) || is.null(true_params$margin2_params)) {
    warning("True params/DGP info/margin params missing: ", basename(sim_data_file), call. = FALSE)
    files_failed <- files_failed + 1
    failure_logs[[length(failure_logs) + 1]] <- tibble(
      file = base_name,
      condition_id = cond_id,
      rep_i = rep_i,
      reason = "true params missing"
    )
    next
  }

  summary_df <- tryCatch(
    {
      summary(fit)$summary %>%
        as.data.frame() %>%
        rownames_to_column("parameter") %>%
        as_tibble() %>%
        select(parameter, post_mean = mean, post_median = `50%`, post_sd = sd, ci_low = `2.5%`, ci_high = `97.5%`, n_eff, Rhat)
    },
    error = function(e) NULL
  )
  if (is.null(summary_df)) {
    files_failed <- files_failed + 1
    failure_logs[[length(failure_logs) + 1]] <- tibble(
      file = base_name,
      condition_id = cond_id,
      rep_i = rep_i,
      reason = "summary extraction failed"
    )
    next
  }

  summary_df <- summary_df %>%
    mutate(
      condition_id = cond_id, rep_i = rep_i, fitted_model_code = fitted_model_code,
      true_value = map_dbl(parameter, ~ get_true_value_var1(.x, true_params, fitted_model_code)),
      bias = post_mean - true_value,
      coverage = ifelse(!is.na(true_value), (true_value >= ci_low & true_value <= ci_high), NA),
      ci_width = ci_high - ci_low,
      rel_bias = ifelse(!is.na(bias) & !is.na(true_value) & abs(true_value) > 1e-6, bias / abs(true_value), NA_real_),
      param_category = map_chr(parameter, get_param_category_var1)
    )
  param_results_list[[length(param_results_list) + 1]] <- summary_df

  # Sampler Diagnostics
  sampler_params <- tryCatch(
    {
      rstan::get_sampler_params(fit, inc_warmup = FALSE)
    },
    error = function(e) NULL
  )
  if (!is.null(sampler_params)) {
    sampler_params_df <- map_dfr(seq_along(sampler_params), ~ as_tibble(sampler_params[[.x]]) %>% mutate(chain = .x))
    total_divergences <- sum(sampler_params_df$divergent__, na.rm = TRUE)
    max_td_setting <- 10
    try(
      {
        max_td_setting <- attr(fit@sim$samples[[1]], "sampler_params")$max_treedepth
      },
      silent = TRUE
    )
    maxdepth_exceeded <- sum(sampler_params_df$treedepth__ >= max_td_setting, na.rm = TRUE)
    avg_eFMI <- mean(map_dbl(sampler_params, ~ calc_e_fmi(.x[, "energy__"])), na.rm = TRUE)

    # Per-chain averages for acceptance statistic and stepsize
    chain_accept <- map_dbl(sampler_params, ~ mean(.x[, "accept_stat__"], na.rm = TRUE))
    chain_steps <- map_dbl(sampler_params, ~ mean(.x[, "stepsize__"], na.rm = TRUE))
    avg_accept_stat <- mean(chain_accept, na.rm = TRUE)
    avg_stepsize <- mean(chain_steps, na.rm = TRUE)

    sampler_results_list[[length(sampler_results_list) + 1]] <- tibble(
      condition_id = cond_id, rep_i = rep_i, fitted_model_code = fitted_model_code,
      divergences = total_divergences, maxdepth_exceeded = maxdepth_exceeded, eFMI = avg_eFMI,
      avg_accept_stat = avg_accept_stat, avg_stepsize = avg_stepsize
    )
  } else {
    sampler_results_list[[length(sampler_results_list) + 1]] <- tibble(
      condition_id = cond_id, rep_i = rep_i, fitted_model_code = fitted_model_code,
      divergences = NA_integer_, maxdepth_exceeded = NA_integer_, eFMI = NA_real_,
      avg_accept_stat = NA_real_, avg_stepsize = NA_real_
    )
  }
  files_processed <- files_processed + 1
}
cat("\n") # Newline after processing loop

cat("Processed", files_processed, "fit files.", files_failed, "failed.\n")

# --- Aggregate & Save Results ---
if (length(param_results_list) > 0) {
  all_param_results <- bind_rows(param_results_list)
  sim_factors <- sim_conditions_df %>% select(condition_id, dgp_copula_type, dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22)
  all_param_results_merged <- all_param_results %>% left_join(sim_factors, by = "condition_id")
  param_output_file <- file.path(RESULTS_DIR, "parameter_summary.rds") # Simplified name
  saveRDS(all_param_results_merged, file = param_output_file)
  cat("Saved parameter results to:", param_output_file, "\n")
} else {
  warning("No parameter results processed.")
}

if (length(sampler_results_list) > 0) {
  all_sampler_results <- bind_rows(sampler_results_list)
  sim_factors <- sim_conditions_df %>% select(condition_id, dgp_copula_type, dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22)
  all_sampler_results_merged <- all_sampler_results %>% left_join(sim_factors, by = "condition_id")
  sampler_output_file <- file.path(RESULTS_DIR, "sampler_summary.rds") # Simplified name
  saveRDS(all_sampler_results_merged, file = sampler_output_file)
  cat("Saved sampler information to:", sampler_output_file, "\n")
} else {
  warning("No sampler results processed.")
}

# Save details about failed fits for debugging
if (length(failure_logs) > 0) {
  failed_df <- bind_rows(failure_logs) %>% mutate(
    file_size = file.info(file.path(FITS_DIR, file))$size
  )
  failure_report <- file.path(RESULTS_DIR, "failed_fit_report.csv")
  write.csv(failed_df, failure_report, row.names = FALSE)
  cat("Saved failed fit report to:", failure_report, "\n")
}

cat("--- Fit Processing Finished ---\n")
