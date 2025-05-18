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

# --- Configuration ---
RESULTS_DIR <- getwd()
DATA_DIR <- file.path(RESULTS_DIR, "../data")
PLOTS_DIR <- file.path(RESULTS_DIR, "plots_convergence")
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
cat("Created directory for plots:", PLOTS_DIR, "\n")
data_dir_abs <- normalizePath(DATA_DIR, mustWork = FALSE)
if (!dir.exists(data_dir_abs)) stop("Data directory not found: ", data_dir_abs)

# --- Options ---
`%||%` <- function(a, b) if (!is.null(a)) a else b
theme_set(theme_bw(base_size = 11))
cat("Setup complete.\n")

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

  sampler_summary_agg <- sampler_joined_df %>%
    filter(!is.na(fitted_model_code)) %>%
    group_by(T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, fit_type, fitted_model_code) %>%
    summarize(
      n_total_fits = n(), n_reps_diverged = safe_summarize(divergences, function(x) sum(x > 0, na.rm = TRUE)),
      pct_reps_diverged = safe_summarize(divergences > 0, mean, na.rm = TRUE) * 100,
      avg_divergences = safe_summarize(divergences, mean, na.rm = TRUE), max_divergences = safe_summarize(divergences, max, na.rm = TRUE),
      n_reps_maxdepth = safe_summarize(maxdepth_exceeded, function(x) sum(x > 0, na.rm = TRUE)),
      pct_reps_maxdepth = safe_summarize(maxdepth_exceeded > 0, mean, na.rm = TRUE) * 100,
      avg_eFMI = safe_summarize(eFMI, mean, na.rm = TRUE), min_eFMI = safe_summarize(eFMI, min, na.rm = TRUE),
      n_reps_low_bfmi = safe_summarize(eFMI, function(x) sum(x < 0.3, na.rm = TRUE)), .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ replace(., !is.finite(.), NA)))

  cat("\n--- Overall Sampler Diagnostics Summary (Aggregated) ---\n")
  print(kable(sampler_summary_agg, digits = 2, row.names = FALSE) %>% kable_styling(font_size = 8))
  cat("-------------------------------------------------------\n")
} else {
  cat("No sampler information available.\n")
}

# --- Generate and Save Aggregated Diagnostic Plots ---
cat("\nGenerating aggregated diagnostic plots...\n")
fill_scale_viridis <- scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fitted Model")
plot_saved_count <- 0

if (nrow(sampler_joined_df) > 0) {
  p_div_agg <- ggplot(sampler_joined_df, aes(x = fit_type, y = divergences, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    labs(title = "Divergences (Aggregated)", x = "Fit Type", y = "Count") +
    theme(legend.position = "bottom")
  ggsave(file.path(PLOTS_DIR, "diag_agg_divergences.png"), p_div_agg, width = 7, height = 5, dpi = 150)
  plot_saved_count <- plot_saved_count + 1

  p_efmi_agg <- ggplot(sampler_joined_df, aes(x = fit_type, y = eFMI, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    geom_hline(yintercept = 0.3, lty = "dashed", color = "red") +
    fill_scale_viridis +
    guides(fill = "none") +
    labs(title = "E-FMI (Aggregated)", x = "Fit Type", y = "E-FMI")
  ggsave(file.path(PLOTS_DIR, "diag_agg_eFMI.png"), p_efmi_agg, width = 7, height = 5, dpi = 150)
  plot_saved_count <- plot_saved_count + 1
}
if (nrow(results_joined_df) > 0) {
  p_rhat_agg <- ggplot(results_joined_df, aes(x = fit_type, y = Rhat, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    geom_hline(yintercept = 1.05, lty = "dashed", color = "red") +
    fill_scale_viridis +
    guides(fill = "none") +
    labs(title = "Rhat (Aggregated)", x = "Fit Type", y = "Rhat") +
    coord_cartesian(ylim = c(0.98, 1.1))
  ggsave(file.path(PLOTS_DIR, "diag_agg_Rhat.png"), p_rhat_agg, width = 7, height = 5, dpi = 150)
  plot_saved_count <- plot_saved_count + 1

  p_neff_agg <- ggplot(results_joined_df, aes(x = fit_type, y = n_eff, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    scale_y_log10(limits = c(1, NA), oob = scales::squish, breaks = scales::log_breaks(n = 5)) +
    annotation_logticks(sides = "l", short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm")) +
    labs(title = "N_eff (Aggregated)", x = "Fit Type", y = "N_eff (log scale)")
  ggsave(file.path(PLOTS_DIR, "diag_agg_Neff.png"), p_neff_agg, width = 7, height = 5, dpi = 150)
  plot_saved_count <- plot_saved_count + 1
}
cat("Saved", plot_saved_count, "aggregated diagnostic plots.\n")

# --- Generate and Save Faceted Diagnostic Plots ---
cat("\nGenerating faceted diagnostic plots...\n")
plot_saved_count <- 0
# Define faceting structure (T vs Alphas + Copula) - adjust if too busy
facet_formula <- T_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac
cat("Using facet formula:", deparse(facet_formula), "\n")

if (nrow(sampler_joined_df) > 0) {
  p_div_facet <- ggplot(sampler_joined_df, aes(x = fit_type, y = divergences, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "Divergences by Condition", x = NULL, y = "Count") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  ggsave(file.path(PLOTS_DIR, "diag_facet_divergences.png"), p_div_facet, width = 12, height = 7, dpi = 150)
  plot_saved_count <- plot_saved_count + 1

  p_efmi_facet <- ggplot(sampler_joined_df, aes(x = fit_type, y = eFMI, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    geom_hline(yintercept = 0.3, lty = "dashed", color = "red") +
    fill_scale_viridis +
    guides(fill = "none") +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "E-FMI by Condition", x = NULL, y = "E-FMI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  ggsave(file.path(PLOTS_DIR, "diag_facet_eFMI.png"), p_efmi_facet, width = 12, height = 7, dpi = 150)
  plot_saved_count <- plot_saved_count + 1
}
if (nrow(results_joined_df) > 0) {
  avg_rhat_per_fit <- results_joined_df %>%
    filter(!is.na(T_fac) & !is.na(dgp_alpha1_fac) & !is.na(dgp_alpha2_fac) & !is.na(dgp_copula_fac)) %>%
    group_by(condition_id, rep_i, fitted_model_code, fit_type, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac) %>%
    summarise(avg_Rhat = mean(Rhat, na.rm = TRUE), .groups = "drop")
  p_rhat_facet <- ggplot(avg_rhat_per_fit, aes(x = fit_type, y = avg_Rhat, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    geom_hline(yintercept = 1.05, lty = "dashed", color = "red") +
    fill_scale_viridis +
    guides(fill = "none") +
    facet_grid(facet_formula, labeller = label_value) +
    coord_cartesian(ylim = c(0.99, 1.06)) +
    labs(title = "Average Rhat by Condition", x = NULL, y = "Average Rhat") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  ggsave(file.path(PLOTS_DIR, "diag_facet_Rhat.png"), p_rhat_facet, width = 12, height = 7, dpi = 150)
  plot_saved_count <- plot_saved_count + 1

  med_neff_per_fit <- results_joined_df %>%
    filter(!is.na(T_fac) & !is.na(dgp_alpha1_fac) & !is.na(dgp_alpha2_fac) & !is.na(dgp_copula_fac)) %>%
    group_by(condition_id, rep_i, fitted_model_code, fit_type, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac) %>%
    summarise(med_n_eff = median(n_eff, na.rm = TRUE), .groups = "drop")
  p_neff_facet <- ggplot(med_neff_per_fit, aes(x = fit_type, y = med_n_eff, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    scale_y_log10(limits = c(10, NA), oob = scales::squish, breaks = scales::log_breaks(n = 4)) +
    annotation_logticks(sides = "l", short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm")) +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "Median N_eff by Condition", x = NULL, y = "Median N_eff (log scale)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  ggsave(file.path(PLOTS_DIR, "diag_facet_Neff.png"), p_neff_facet, width = 12, height = 7, dpi = 150)
  plot_saved_count <- plot_saved_count + 1
}
cat("Saved", plot_saved_count, "faceted diagnostic plots.\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\n--- Convergence Analysis Script Finished ---\n")
