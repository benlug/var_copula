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
  # violin and box plots across all conditions
  p_div_agg <- ggplot(sampler_joined_df, aes(x = fit_type, y = divergences, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    labs(title = "Divergences (Aggregated)", x = "Fit Type", y = "Count") +
    theme(legend.position = "bottom")
  print(p_div_agg)

  p_efmi_agg <- ggplot(sampler_joined_df, aes(x = fit_type, y = eFMI, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    geom_hline(yintercept = 0.3, lty = "dashed", color = "red") +
    fill_scale_viridis +
    guides(fill = "none") +
    labs(title = "E-FMI (Aggregated)", x = "Fit Type", y = "E-FMI")
  print(p_efmi_agg)

  p_accept_agg <- ggplot(sampler_joined_df, aes(x = fit_type, y = avg_accept_stat, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    labs(title = "Avg Acceptance Stat (Aggregated)", x = "Fit Type", y = "Acceptance")
  print(p_accept_agg)

  p_stepsize_agg <- ggplot(sampler_joined_df, aes(x = fit_type, y = avg_stepsize, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    labs(title = "Avg Stepsize (Aggregated)", x = "Fit Type", y = "Stepsize")
  print(p_stepsize_agg)
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
  print(p_rhat_agg)

  p_neff_agg <- ggplot(results_joined_df, aes(x = fit_type, y = n_eff, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    scale_y_log10(limits = c(1, NA), oob = scales::squish, breaks = scales::log_breaks(n = 5)) +
    annotation_logticks(sides = "l", short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm")) +
    labs(title = "N_eff (Aggregated)", x = "Fit Type", y = "N_eff (log scale)")
  print(p_neff_agg)
}
cat("Aggregated diagnostic plots generated.\n")

# --- Generate and Save Faceted Diagnostic Plots ---
cat("\nGenerating faceted diagnostic plots...\n")
# Define faceting structure including Tau
facet_formula <- T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac
cat("Using facet formula:", deparse(facet_formula), "\n")

if (nrow(sampler_joined_df) > 0) {
  p_div_facet <- ggplot(sampler_joined_df, aes(x = fit_type, y = divergences, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "Divergences by Condition", x = NULL, y = "Count") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  print(p_div_facet)

  p_efmi_facet <- ggplot(sampler_joined_df, aes(x = fit_type, y = eFMI, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    geom_hline(yintercept = 0.3, lty = "dashed", color = "red") +
    fill_scale_viridis +
    guides(fill = "none") +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "E-FMI by Condition", x = NULL, y = "E-FMI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  print(p_efmi_facet)

  p_accept_facet <- ggplot(sampler_joined_df, aes(x = fit_type, y = avg_accept_stat, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "Avg Acceptance Stat by Condition", x = NULL, y = "Acceptance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  print(p_accept_facet)

  p_stepsize_facet <- ggplot(sampler_joined_df, aes(x = fit_type, y = avg_stepsize, fill = fitted_model_code)) +
    geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
    fill_scale_viridis +
    guides(fill = "none") +
    facet_grid(facet_formula, labeller = label_value) +
    labs(title = "Avg Stepsize by Condition", x = NULL, y = "Stepsize") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(size = 7))
  print(p_stepsize_facet)
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
  print(p_rhat_facet)

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
  print(p_neff_facet)
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
    p_div <- ggplot(cond_samp, aes(x = fit_type, y = divergences, fill = fitted_model_code)) +
      geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
      fill_scale_viridis +
      labs(title = "Divergences", x = "Fit Type", y = "Count") +
      theme(legend.position = "bottom")
    print(p_div)

    p_efmi <- ggplot(cond_samp, aes(x = fit_type, y = eFMI, fill = fitted_model_code)) +
      geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
      geom_hline(yintercept = 0.3, lty = "dashed", color = "red") +
      fill_scale_viridis +
      guides(fill = "none") +
      labs(title = "E-FMI", x = "Fit Type", y = "E-FMI")
    print(p_efmi)

    p_accept <- ggplot(cond_samp, aes(x = fit_type, y = avg_accept_stat, fill = fitted_model_code)) +
      geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
      fill_scale_viridis +
      guides(fill = "none") +
      labs(title = "Avg Acceptance", x = "Fit Type", y = "Acceptance")
    print(p_accept)

    p_step <- ggplot(cond_samp, aes(x = fit_type, y = avg_stepsize, fill = fitted_model_code)) +
      geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
      fill_scale_viridis +
      guides(fill = "none") +
      labs(title = "Avg Stepsize", x = "Fit Type", y = "Stepsize")
    print(p_step)
  }
  if (nrow(cond_res) > 0) {
    avg_rhat <- cond_res %>%
      group_by(condition_id, rep_i, fitted_model_code, fit_type) %>%
      summarise(avg_Rhat = mean(Rhat, na.rm = TRUE), .groups = "drop")
    p_rhat <- ggplot(avg_rhat, aes(x = fit_type, y = avg_Rhat, fill = fitted_model_code)) +
      geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
      geom_hline(yintercept = 1.05, lty = "dashed", color = "red") +
      fill_scale_viridis +
      guides(fill = "none") +
      coord_cartesian(ylim = c(0.99, 1.06)) +
      labs(title = "Average Rhat", x = "Fit Type", y = "Avg Rhat")
    print(p_rhat)

    med_neff <- cond_res %>%
      group_by(condition_id, rep_i, fitted_model_code, fit_type) %>%
      summarise(med_n_eff = median(n_eff, na.rm = TRUE), .groups = "drop")
    p_neff <- ggplot(med_neff, aes(x = fit_type, y = med_n_eff, fill = fitted_model_code)) +
      geom_violin(trim = F, alpha = 0.6, na.rm = T, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.8, na.rm = T) +
      fill_scale_viridis +
      guides(fill = "none") +
      scale_y_log10(limits = c(10, NA), oob = scales::squish, breaks = scales::log_breaks(n = 4)) +
      annotation_logticks(sides = "l") +
      labs(title = "Median N_eff", x = "Fit Type", y = "Median N_eff (log)")
    print(p_neff)
  }

  # Pairs plot for the replicate with most divergences
  div_run <- cond_samp %>%
    filter(divergences > 0) %>%
    arrange(desc(divergences)) %>%
    slice(1)
  if (nrow(div_run) == 1) {
    fit_path <- file.path(FITS_DIR, sprintf("fit_%s_cond%03d_rep%03d.rds", div_run$fitted_model_code, condition_id, div_run$rep_i))
    fit_obj <- safe_read_stanfit(fit_path)
    if (!is.null(fit_obj)) {
      # bayesplot::mcmc_pairs expects a NUTS parameter data frame with
      # columns Chain, Iteration, Parameter, Value. Construct this from the
      # sampler parameters, focusing on divergences.
      sampler_list <- rstan::get_sampler_params(fit_obj, inc_warmup = FALSE)
      np_df <- purrr::imap_dfr(sampler_list, ~ {
        tibble(
          Chain = .y,
          Iteration = seq_len(nrow(.x)),
          Parameter = "divergent__",
          Value = .x[, "divergent__"]
        )
      })
      arr <- rstan::extract(fit_obj, permuted = FALSE)
      par_names <- dimnames(arr)$parameters
      sel_pars <- head(par_names[!grepl("^lp__", par_names)], 5)
      if (length(sel_pars) >= 2) {
        p_pairs <- bayesplot::mcmc_pairs(arr, pars = sel_pars, np = np_df)
        print(p_pairs)
      }
    }
  }
  grDevices::dev.off()
  cat("Saved", pdf_file, "\n")
})

cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\n--- Convergence Analysis Script Finished ---\n")
