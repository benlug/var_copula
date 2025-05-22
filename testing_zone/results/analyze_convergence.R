###########################################################################
# analyze_convergence.R
#
# Analyzes MCMC convergence diagnostics from VAR(1) Copula simulation.
# Handles separate alpha1/alpha2. Simplified naming.
# Assumes run from 'results' directory.
###########################################################################

# --- Libraries ---
library(rstan)
library(tidyverse)
library(knitr)
library(kableExtra)
library(forcats)
library(ggplot2)
library(this.path)
library(purrr)
library(stringr)
library(bayesplot)
library(grid)

# --- Configuration ---
RESULTS_DIR <- getwd()
DATA_DIR <- file.path(RESULTS_DIR, "../data")
PLOTS_DIR <- file.path(RESULTS_DIR, "plots_convergence")
FITS_DIR <- file.path(RESULTS_DIR, "../fits")
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
cat("Created directory for plots:", PLOTS_DIR, "\n")
data_dir_abs <- normalizePath(DATA_DIR, mustWork = FALSE)
if (!dir.exists(data_dir_abs)) stop("Data directory not found: ", data_dir_abs)

# --- Options ---
`%||%` <- function(a, b) if (!is.null(a)) a else b
theme_set(theme_bw(base_size = 11))
cat("Setup complete.\n")

safe_read_stanfit <- function(filename) {
  if (!file.exists(filename)) {
    return(NULL)
  }
  tryCatch(
    {
      fit <- readRDS(filename)
      if (inherits(fit, "stanfit")) fit else NULL
    },
    error = function(e) NULL
  )
}

# --- Load Simulation Conditions ---
sim_conds_file <- file.path(data_dir_abs, "sim_conditions.rds")
if (!file.exists(sim_conds_file)) stop("Sim conditions file not found: ", sim_conds_file)
sim_conditions_df <- readRDS(sim_conds_file) %>%
  mutate(condition_id = as.integer(condition_id)) %>%
  select(condition_id, dgp_copula_type, dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22, correct_fit_model_code) %>%
  distinct(condition_id, .keep_all = TRUE)
cat("Loaded simulation conditions.\n")

# --- Load Parameter Results (for Rhat/Neff) ---
param_results_file <- file.path(RESULTS_DIR, "parameter_summary.rds")
if (!file.exists(param_results_file)) stop("Parameter results file not found: ", param_results_file)
results_df <- readRDS(param_results_file) %>% mutate(condition_id = as.integer(condition_id))
cat("Loaded parameter results.\n")

# --- Load Sampler Info ---
sampler_info_file <- file.path(RESULTS_DIR, "sampler_summary.rds")
# Corrected if/else structure for loading sampler info
if (!file.exists(sampler_info_file)) {
  warning("Sampler info file not found.")
  sampler_info_df <- data.frame()
} else {
  sampler_info_df <- readRDS(sampler_info_file) %>% mutate(condition_id = as.integer(condition_id))
}
cat("Loaded sampler info.\n")

# --- Basic Checks ---
if (nrow(results_df) == 0 && nrow(sampler_info_df) == 0) stop("No results or sampler info found.")
cat("Results or sampler info loaded.\n")

# --- Prepare Factors ---
fitted_model_levels <- c("NG", "NC", "SG", "SC")
dgp_copula_levels <- sort(unique(c(results_df$dgp_copula_type, sampler_info_df$dgp_copula_type)))
dgp_alpha1_levels <- sort(unique(c(results_df$dgp_alpha1, sampler_info_df$dgp_alpha1)))
dgp_alpha2_levels <- sort(unique(c(results_df$dgp_alpha2, sampler_info_df$dgp_alpha2)))
dgp_tau_levels <- sort(unique(c(results_df$dgp_tau, sampler_info_df$dgp_tau)))
T_levels <- sort(unique(c(results_df$T, sampler_info_df$T)))
master_param_order <- c("phi11", "phi12", "phi21", "phi22", "mu[1]", "mu[2]", "rho", "theta", "tau", "sigma[1]", "sigma[2]", "xi[1]", "xi[2]", "omega[1]", "omega[2]", "alpha[1]", "alpha[2]")

# --- Factorize Parameter Results ---
if (nrow(results_df) > 0) {
  param_levels <- intersect(master_param_order, unique(as.character(results_df$parameter)))
  results_joined_df <- results_df %>%
    left_join(select(sim_conditions_df, condition_id, correct_fit_model_code), by = "condition_id") %>%
    mutate(
      dgp_model_code = correct_fit_model_code,
      fit_type = ifelse(fitted_model_code == dgp_model_code, "Correct", "Standard")
    ) %>%
    mutate(
      parameter = factor(parameter, levels = param_levels),
      fitted_model_code = factor(fitted_model_code, levels = fitted_model_levels),
      fit_type = factor(fit_type, levels = c("Correct", "Standard")),
      T_fac = factor(paste0("T=", T), levels = paste0("T=", T_levels)),
      dgp_alpha1_fac = factor(paste0("A1=", dgp_alpha1), levels = paste0("A1=", dgp_alpha1_levels)),
      dgp_alpha2_fac = factor(paste0("A2=", dgp_alpha2), levels = paste0("A2=", dgp_alpha2_levels)),
      dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)), levels = paste0("Cop=", str_to_title(dgp_copula_levels))),
      dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
    ) %>%
    mutate(parameter = fct_drop(parameter))
  param_levels <- levels(results_joined_df$parameter)
} else {
  results_joined_df <- data.frame()
  param_levels <- character(0)
}

# --- Factorize Sampler Info ---
if (nrow(sampler_info_df) > 0) {
  sampler_joined_df <- sampler_info_df %>%
    left_join(select(sim_conditions_df, condition_id, correct_fit_model_code), by = "condition_id") %>%
    mutate(
      dgp_model_code = correct_fit_model_code,
      fit_type = ifelse(fitted_model_code == dgp_model_code, "Correct", "Standard")
    ) %>%
    mutate(
      fitted_model_code = factor(fitted_model_code, levels = fitted_model_levels),
      fit_type = factor(fit_type, levels = c("Correct", "Standard")),
      T_fac = factor(paste0("T=", T), levels = paste0("T=", T_levels)),
      dgp_alpha1_fac = factor(paste0("A1=", dgp_alpha1), levels = paste0("A1=", dgp_alpha1_levels)),
      dgp_alpha2_fac = factor(paste0("A2=", dgp_alpha2), levels = paste0("A2=", dgp_alpha2_levels)),
      dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)), levels = paste0("Cop=", str_to_title(dgp_copula_levels))),
      dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
    )
} else {
  sampler_joined_df <- data.frame()
}
cat("Data preparation complete.\n")

# --- Calculate Overall Sampler Diagnostics Summary ---
cat("\nCalculating overall sampler diagnostics summary...\n")
if (nrow(sampler_joined_df) > 0) {
  safe_summarize <- function(x, fun, ...) { # Simplified
    x_finite <- x[is.finite(x)]
    if (length(x_finite) == 0) {
      return(NA)
    }
    tryCatch(fun(x_finite, ...), error = function(e) NA)
  }

  # summarize sampler diagnostics across replications
  sampler_summary_agg <- sampler_joined_df %>%
    filter(!is.na(fitted_model_code)) %>%
    group_by(condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, fit_type, fitted_model_code) %>%
    summarize(
      n_total_fits = n(), n_reps_diverged = safe_summarize(divergences, function(x) sum(x > 0, na.rm = TRUE)),
      pct_reps_diverged = safe_summarize(divergences > 0, mean, na.rm = TRUE) * 100,
      avg_divergences = safe_summarize(divergences, mean, na.rm = TRUE), max_divergences = safe_summarize(divergences, max, na.rm = TRUE),
      n_reps_maxdepth = safe_summarize(maxdepth_exceeded, function(x) sum(x > 0, na.rm = TRUE)),
      pct_reps_maxdepth = safe_summarize(maxdepth_exceeded > 0, mean, na.rm = TRUE) * 100,
      avg_eFMI = safe_summarize(eFMI, mean, na.rm = TRUE), min_eFMI = safe_summarize(eFMI, min, na.rm = TRUE),
      n_reps_low_bfmi = safe_summarize(eFMI, function(x) sum(x < 0.3, na.rm = TRUE)),
      avg_accept_stat = safe_summarize(avg_accept_stat, mean, na.rm = TRUE),
      avg_stepsize = safe_summarize(avg_stepsize, mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ replace(., !is.finite(.), NA)))

  cat("\n--- Overall Sampler Diagnostics Summary (Aggregated) ---\n")
  print(kable(sampler_summary_agg, digits = 2, row.names = FALSE) %>% kable_styling(font_size = 8))
  cat("-------------------------------------------------------\n")

  # Save aggregated summary to CSV
  summary_csv <- file.path(RESULTS_DIR, "convergence_summary.csv")
  write.csv(sampler_summary_agg, summary_csv, row.names = FALSE)
  cat("Saved convergence summary to:", summary_csv, "\n")
} else {
  cat("No sampler information available.\n")
}

# --- Generate and Save Aggregated Diagnostic Plots ---
cat("\nGenerating aggregated diagnostic plots...\n")
fill_scale_viridis <- scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fitted Model")

if (nrow(sampler_joined_df) > 0) {
  # histogram of divergences across all conditions
  p_div_agg <- ggplot(sampler_joined_df, aes(x = divergences, fill = fitted_model_code)) +
    geom_histogram(binwidth = 1, colour = "black", alpha = 0.7, position = "dodge") +
    fill_scale_viridis +
    labs(title = "Histogram of Divergent Transitions", x = "Number of divergences", y = "Replications") +
    theme(legend.position = "bottom")
  print(p_div_agg)
}
if (nrow(results_joined_df) > 0) {
  # no additional aggregated plots required
}
cat("Aggregated diagnostic plots generated.\n")

# --- Generate and Save Faceted Diagnostic Plots ---
cat("\nGenerating faceted diagnostic plots...\n")
# Define faceting structure including Tau
facet_formula <- T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac
cat("Using facet formula:", deparse(facet_formula), "\n")

if (nrow(sampler_joined_df) > 0) {
  p_div_facet <- ggplot(sampler_joined_df, aes(x = divergences, fill = fitted_model_code)) +
    geom_histogram(binwidth = 1, colour = "black", position = "dodge") +
    fill_scale_viridis +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "Divergent Transitions by Condition", x = "Number of divergences", y = "Replications") +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 8), strip.text = element_text(size = 7))
  print(p_div_facet)
}
if (nrow(results_joined_df) > 0) {
  # no additional faceted plots required
}
cat("Faceted diagnostic plots generated.\n")

# --- Per-Condition PDF Plots ---
condition_info_df <- sim_conditions_df %>%
  mutate(
    T_fac = factor(paste0("T=", T), levels = paste0("T=", T_levels)),
    dgp_alpha1_fac = factor(paste0("A1=", dgp_alpha1), levels = paste0("A1=", dgp_alpha1_levels)),
    dgp_alpha2_fac = factor(paste0("A2=", dgp_alpha2), levels = paste0("A2=", dgp_alpha2_levels)),
    dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)), levels = paste0("Cop=", str_to_title(dgp_copula_levels))),
    dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
  ) %>%
  select(condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, phi11, phi12, phi21, phi22)

cat("Generating per-condition convergence PDFs...\n")
purrr::pwalk(condition_info_df, function(condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, phi11, phi12, phi21, phi22) {
  cond_res <- results_joined_df %>% filter(condition_id == !!condition_id)
  cond_samp <- sampler_joined_df %>% filter(condition_id == !!condition_id)
  if (nrow(cond_res) == 0 && nrow(cond_samp) == 0) {
    return()
  }
  pdf_file <- file.path(PLOTS_DIR, sprintf("convergence_cond_%03d.pdf", condition_id))
  grDevices::pdf(pdf_file, width = 7, height = 5)
  grid::grid.newpage()
  grid::grid.text(sprintf("Condition %03d\n%s %s %s %s %s\nphi11=%.2f phi12=%.2f phi21=%.2f phi22=%.2f", condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, phi11, phi12, phi21, phi22), gp = grid::gpar(cex = 0.9))
  if (nrow(cond_samp) > 0) {
    p_hist <- ggplot(cond_samp, aes(x = divergences, fill = fitted_model_code)) +
      geom_histogram(binwidth = 1, colour = "black", position = "dodge") +
      fill_scale_viridis +
      labs(title = "Divergent Transitions", x = "Number of divergences", y = "Replications") +
      theme(legend.position = "bottom")
    print(p_hist)
  }
  grDevices::dev.off()
  cat("Saved", pdf_file, "\n")
})

cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\n--- Convergence Analysis Script Finished ---\n")
