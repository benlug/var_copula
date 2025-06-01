###########################################################################
# analyze_parameter_recovery.R
#
# Analyzes parameter recovery from VAR(1) Copula simulation (separate alphas).
# Creates multiple PDFs:
#   1) Main metrics (rel_bias, post_sd)
#   2) Detailed metrics (raw bias, RMSE)
#   3) Special analysis: avg_bias, coverage, sd_est (all reps vs no divergences),
#      plus an alpha-parameters scatter plot of avg_est vs. avg_rel_bias
#
# Assumes run from 'results' directory.
###########################################################################

# --- Libraries ---
library(tidyverse)
library(knitr)
library(kableExtra)
library(forcats)
library(ggplot2)
library(this.path)
library(purrr)
library(stringr)
library(grid)

# --- Configuration ---
RESULTS_DIR <- tryCatch(this.dir(), error = function(e) getwd())
DATA_DIR <- file.path(RESULTS_DIR, "../data")
PLOTS_DIR <- file.path(RESULTS_DIR, "plots_param_recovery")
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
cat("Created directory for plots:", PLOTS_DIR, "\n")

data_dir_abs <- normalizePath(DATA_DIR, mustWork = FALSE)
if (!dir.exists(data_dir_abs)) stop("Data directory not found: ", data_dir_abs)

# --- Options ---
`%||%` <- function(a, b) if (!is.null(a)) a else b
theme_set(theme_bw(base_size = 10))
cat("Setup complete.\n")

# --- Load Data ---
sim_conds_file <- file.path(data_dir_abs, "sim_conditions.rds")
if (!file.exists(sim_conds_file)) stop("Sim conditions file not found: ", sim_conds_file)
sim_conditions_df <- readRDS(sim_conds_file) %>%
  mutate(condition_id = as.integer(condition_id)) %>%
  select(condition_id, dgp_copula_type, dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22, correct_fit_model_code) %>%
  distinct(condition_id, .keep_all = TRUE)

param_results_file <- file.path(RESULTS_DIR, "parameter_summary.rds")
if (!file.exists(param_results_file)) stop("Parameter results file not found: ", param_results_file)
results_df <- readRDS(param_results_file) %>% mutate(condition_id = as.integer(condition_id))

sampler_info_file <- file.path(RESULTS_DIR, "sampler_summary.rds")
if (!file.exists(sampler_info_file)) {
  warning("Sampler info file not found.")
  sampler_info_df <- data.frame()
} else {
  sampler_info_df <- readRDS(sampler_info_file) %>% mutate(condition_id = as.integer(condition_id))
}

if (nrow(results_df) == 0) stop("Parameter results data frame is empty.")
cat("Data loaded successfully.\n")

# --- Define Factor Levels and Order ---
master_param_order <- c(
  "phi11", "phi12", "phi21", "phi22",
  "mu[1]", "mu[2]", "rho", "theta", "tau",
  "sigma[1]", "sigma[2]", "xi[1]", "xi[2]",
  "omega[1]", "omega[2]", "alpha[1]", "alpha[2]"
)
category_order <- c("VAR Coeffs", "Intercepts", "Copula Params", "Residual Params", "Other")
fitted_model_levels <- c("NG", "NC", "SG", "SC")

# --- Determine DGP factor levels DYNAMICALLY ---
dgp_copula_levels <- sort(unique(results_df$dgp_copula_type))
dgp_alpha1_levels <- sort(unique(results_df$dgp_alpha1))
dgp_alpha2_levels <- sort(unique(results_df$dgp_alpha2))
dgp_tau_levels <- sort(unique(results_df$dgp_tau))
T_levels <- sort(unique(results_df$T))
cat("Determined factor levels from data.\n")

# --- Prepare Data for Plotting ---
# Make sure we have rep_i in results_df
if (!"rep_i" %in% names(results_df)) {
  stop("results_df must contain a 'rep_i' column (replicate index).")
}

actual_params_in_data <- unique(as.character(results_df$parameter))
param_levels <- intersect(master_param_order, actual_params_in_data)
actual_categories_in_data <- unique(as.character(results_df$param_category))
category_order <- intersect(category_order, actual_categories_in_data)

results_joined_df <- results_df %>%
  filter(!is.na(parameter)) %>%
  mutate(
    dgp_model_code = paste0(
      toupper(substr(
        ifelse(abs(dgp_alpha1) < 1e-6 & abs(dgp_alpha2) < 1e-6, "normal", "skewnormal"),
        1, 1
      )),
      toupper(substr(dgp_copula_type, 1, 1))
    ),
    fit_type = ifelse(fitted_model_code == dgp_model_code, "Correct", "Standard")
  ) %>%
  mutate(
    parameter = factor(parameter, levels = param_levels),
    param_category = factor(replace_na(param_category, "Other"), levels = category_order),
    fitted_model_code = factor(fitted_model_code, levels = fitted_model_levels),
    fit_type = factor(fit_type, levels = c("Correct", "Standard")),
    T_fac = factor(paste0("T=", T), levels = paste0("T=", T_levels)),
    dgp_alpha1_fac = factor(paste0("A1=", dgp_alpha1), levels = paste0("A1=", dgp_alpha1_levels)),
    dgp_alpha2_fac = factor(paste0("A2=", dgp_alpha2), levels = paste0("A2=", dgp_alpha2_levels)),
    dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)),
      levels = paste0("Cop=", str_to_title(dgp_copula_levels))
    ),
    dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
  )

# --- Factorize Sampler Info ---
if (nrow(sampler_info_df) > 0) {
  if (!"rep_i" %in% names(sampler_info_df)) {
    stop("sampler_info_df must contain a 'rep_i' column to merge on replicate index.")
  }

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
      dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)),
        levels = paste0("Cop=", str_to_title(dgp_copula_levels))
      ),
      dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
    )
} else {
  sampler_joined_df <- data.frame()
}
cat("Data preparation complete.\n")

# --- Helper Plotting Functions ---
plot_boxplot_dist_agg <- function(data, y_var, category_name, title_suffix, y_label,
                                  hline_val = NULL, param_levels_arg, explanation_text = "") {
  # Creates a boxplot over replications for metric y_var, restricted to param_category
  # 'explanation_text' appended to the subtitle for interpretation hints.

  y_sym <- sym(y_var)
  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))

  # If no rows or no param levels, skip
  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
    return(NULL)
  }

  # NEW: Check if all Y values are NA
  valid_y_count <- sum(!is.na(data_plot[[y_var]]))
  if (valid_y_count == 0) {
    # No valid numeric data to plot
    return(NULL)
  }

  p <- ggplot(data_plot, aes(x = parameter, y = !!y_sym, fill = fit_type)) +
    geom_boxplot(
      position = position_dodge(width = 0.85),
      outlier.shape = ".",
      outlier.alpha = 0.5,
      alpha = 0.7,
      width = 0.7,
      na.rm = TRUE
    ) +
    {
      if (!is.null(hline_val)) {
        geom_hline(yintercept = hline_val, color = "red", linetype = "dashed")
      }
    } +
    scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fit Type") +
    scale_x_discrete(limits = param_levels_arg, drop = FALSE) +
    labs(
      title = paste(category_name, "-", title_suffix),
      subtitle = explanation_text,
      x = NULL,
      y = y_label
    ) +
    facet_grid(
      T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac,
      scales = "free_y",
      labeller = label_value
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10),
      strip.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      strip.text = element_text(size = 7)
    )

  return(p)
}

plot_rmse_bar_agg <- function(data, category_name, title_suffix, param_levels_arg, explanation_text = "") {
  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))

  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
    return(NULL)
  }

  # Check if all rmse are NA
  if (all(is.na(data_plot$rmse))) {
    return(NULL)
  }

  p <- ggplot(data_plot, aes(x = parameter, y = rmse, fill = fit_type)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "grey30",
      alpha = 0.8,
      na.rm = TRUE
    ) +
    scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fit Type") +
    scale_x_discrete(limits = param_levels_arg, drop = FALSE) +
    labs(
      title      = paste(category_name, "-", title_suffix),
      subtitle   = explanation_text,
      x          = NULL,
      y          = "RMSE"
    ) +
    facet_grid(
      T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac,
      scales   = "free_y",
      labeller = label_value
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10),
      strip.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      strip.text = element_text(size = 7)
    )

  return(p)
}

plot_metric_bar_agg <- function(data, metric, category_name, title_suffix,
                                metric_label, param_levels_arg, explanation_text = "") {
  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))

  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
    return(NULL)
  }

  # Check if all metric are NA
  if (all(is.na(data_plot[[metric]]))) {
    return(NULL)
  }

  met_sym <- sym(metric)

  p <- ggplot(data_plot, aes(x = parameter, y = !!met_sym, fill = fit_type)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "grey30",
      alpha = 0.8,
      na.rm = TRUE
    ) +
    scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fit Type") +
    scale_x_discrete(limits = param_levels_arg, drop = FALSE) +
    labs(
      title      = paste(category_name, "-", title_suffix),
      subtitle   = explanation_text,
      x          = NULL,
      y          = metric_label
    ) +
    facet_grid(
      T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac,
      scales   = "free_y",
      labeller = label_value
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10),
      strip.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      strip.text = element_text(size = 7)
    )

  return(p)
}

# --- Calculate Summary Statistics Over All Replications ---
cat("Calculating aggregated summary statistics...\n")
summary_df <- results_joined_df %>%
  filter(!is.na(parameter) & !is.na(fitted_model_code) & !is.na(param_category)) %>%
  group_by(
    condition_id, rep_i, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
    fit_type, fitted_model_code, parameter, param_category
  ) %>%
  summarize(
    bias         = first(bias),
    coverage     = first(coverage),
    post_sd      = first(post_sd),
    rel_bias     = first(rel_bias),
    post_mean    = first(post_mean),
    sd_est       = NA_real_,
    .groups      = "drop"
  )

summary_agg_df <- summary_df %>%
  group_by(
    T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
    fit_type, fitted_model_code, parameter, param_category
  ) %>%
  summarize(
    n_reps_param = n(),
    true_value = NA_real_,
    avg_bias = mean(bias, na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    avg_post_sd = mean(post_sd, na.rm = TRUE),
    sd_est = sd(
      results_df$post_mean[results_df$parameter == first(parameter) &
        results_df$condition_id == first(condition_id) &
        results_df$fitted_model_code == first(fitted_model_code)],
      na.rm = TRUE
    ),
    avg_rel_bias = mean(rel_bias, na.rm = TRUE),
    .groups = "drop"
  )

# Merge sampler info for no-divergences aggregator
if (nrow(sampler_joined_df) > 0) {
  sampler_diverge <- sampler_joined_df %>%
    select(condition_id, rep_i, fitted_model_code, divergences) %>%
    distinct()

  summary_df_nodiv <- summary_df %>%
    left_join(sampler_diverge,
      by = c("condition_id", "rep_i", "fitted_model_code")
    ) %>%
    filter(is.na(divergences) | divergences == 0)

  summary_agg_nodiv_df <- summary_df_nodiv %>%
    group_by(
      T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
      fit_type, fitted_model_code, parameter, param_category
    ) %>%
    summarize(
      n_reps_param = n(),
      avg_bias     = mean(bias, na.rm = TRUE),
      coverage     = mean(coverage, na.rm = TRUE),
      avg_post_sd  = mean(post_sd, na.rm = TRUE),
      sd_est       = sd(post_mean, na.rm = TRUE),
      avg_rel_bias = mean(rel_bias, na.rm = TRUE),
      .groups      = "drop"
    )
} else {
  summary_agg_nodiv_df <- NULL
}

write.csv(summary_agg_df, file.path(RESULTS_DIR, "summary_aggregated_allreps.csv"), row.names = FALSE)
if (!is.null(summary_agg_nodiv_df)) {
  write.csv(summary_agg_nodiv_df, file.path(RESULTS_DIR, "summary_aggregated_nodiv.csv"), row.names = FALSE)
}
cat("Saved aggregated summary statistics.\n")

# -----------------------------------------------------------------------------
#  1) MAIN PDF: relative bias + posterior SD
# -----------------------------------------------------------------------------
pdf_main <- file.path(PLOTS_DIR, "parameter_recovery_faceted_main.pdf")
grDevices::pdf(pdf_main, width = 11, height = 6)
cat("Creating MAIN PDF with relative bias + posterior SD...\n")

results_joined_df_main <- results_joined_df %>%
  mutate(
    rel_bias_capped = ifelse(
      !is.na(rel_bias), pmax(-2, pmin(2, rel_bias)), NA_real_
    )
  )

main_plots <- list(
  list(
    metric         = "rel_bias_capped",
    label          = "Relative Bias (capped ±2)",
    hline          = 0,
    explanation    = "Relative bias = (Estimate - True) / |True|. Zero = no bias.",
    type           = "box"
  ),
  list(
    metric         = "post_sd",
    label          = "Posterior SD",
    explanation    = "Posterior SD from each chain's inference. Higher = more uncertainty.",
    type           = "box"
  )
)

for (cat_name in category_order) {
  category_params <- intersect(
    master_param_order,
    unique(as.character(results_joined_df_main %>% filter(param_category == cat_name) %>% pull(parameter)))
  )
  if (length(category_params) == 0) next

  for (plot_info in main_plots) {
    plot_obj <- NULL
    if (plot_info$type == "box") {
      plot_obj <- plot_boxplot_dist_agg(
        data = results_joined_df_main,
        y_var = plot_info$metric,
        category_name = cat_name,
        title_suffix = "All Conditions",
        y_label = plot_info$label,
        hline_val = plot_info$hline,
        param_levels_arg = category_params,
        explanation_text = plot_info$explanation
      )
    }
    if (!is.null(plot_obj)) print(plot_obj)
  }
}
dev.off()
cat("Saved main PDF:", pdf_main, "\n")

# -----------------------------------------------------------------------------
#  2) DETAILED PDF: raw bias + RMSE
# -----------------------------------------------------------------------------
pdf_detail <- file.path(PLOTS_DIR, "parameter_recovery_faceted_detailed.pdf")
grDevices::pdf(pdf_detail, width = 11, height = 6)
cat("Creating DETAILED PDF with raw bias + RMSE...\n")

detailed_plots <- list(
  list(
    metric         = "bias",
    label          = "Raw Bias",
    hline          = 0,
    explanation    = "Bias = (Estimate - True). Zero = no bias.",
    type           = "box"
  )
  # You can add RMSE aggregator here if you want
)

for (cat_name in category_order) {
  category_params <- intersect(
    master_param_order,
    unique(as.character(results_joined_df %>% filter(param_category == cat_name) %>% pull(parameter)))
  )
  if (length(category_params) == 0) next

  for (plot_info in detailed_plots) {
    plot_obj <- NULL
    if (plot_info$type == "box") {
      plot_obj <- plot_boxplot_dist_agg(
        data = results_joined_df,
        y_var = plot_info$metric,
        category_name = cat_name,
        title_suffix = "All Conditions",
        y_label = plot_info$label,
        hline_val = plot_info$hline,
        param_levels_arg = category_params,
        explanation_text = plot_info$explanation
      )
    }
    if (!is.null(plot_obj)) print(plot_obj)
  }
}
dev.off()
cat("Saved detailed PDF:", pdf_detail, "\n")

# -----------------------------------------------------------------------------
#  3) SPECIAL PDF: avg_bias, coverage, sd_est (all reps vs no divergences)
#     + alpha-parameter scatter (avg_est vs avg_rel_bias)
# -----------------------------------------------------------------------------
pdf_special <- file.path(PLOTS_DIR, "parameter_recovery_faceted_special.pdf")
grDevices::pdf(pdf_special, width = 11, height = 6)
cat("Creating SPECIAL PDF with avg_bias, coverage, sd_est, alpha scatter...\n")

if (!is.null(summary_agg_df)) {
  cat("  -> Plotting aggregator metrics for ALL reps...\n")
  special_metrics <- list(
    list(
      metric = "avg_bias", label = "Average Bias",
      explanation = "Mean difference (Estimate - True) over all replications."
    ),
    list(
      metric = "coverage", label = "Coverage",
      explanation = "Proportion of 95% CIs containing the true value. Ideal near 0.95."
    ),
    list(
      metric = "sd_est", label = "sd_est",
      explanation = "Std. deviation of posterior means across replications (empirical spread)."
    )
  )

  for (cat_name in category_order) {
    category_params <- intersect(
      master_param_order,
      unique(as.character(summary_agg_df %>% filter(param_category == cat_name) %>% pull(parameter)))
    )
    if (length(category_params) == 0) next

    for (sm in special_metrics) {
      plt <- plot_metric_bar_agg(
        data = summary_agg_df,
        metric = sm$metric,
        category_name = cat_name,
        title_suffix = "All Reps",
        metric_label = sm$label,
        param_levels_arg = category_params,
        explanation_text = sm$explanation
      )
      if (!is.null(plt)) print(plt)
    }
  }

  if (!is.null(summary_agg_nodiv_df)) {
    cat("  -> Plotting aggregator metrics for NO DIVERGENCES...\n")
    for (cat_name in category_order) {
      category_params <- intersect(
        master_param_order,
        unique(as.character(summary_agg_nodiv_df %>% filter(param_category == cat_name) %>% pull(parameter)))
      )
      if (length(category_params) == 0) next

      for (sm in special_metrics) {
        plt <- plot_metric_bar_agg(
          data = summary_agg_nodiv_df,
          metric = sm$metric,
          category_name = cat_name,
          title_suffix = "No Divergences",
          metric_label = sm$label,
          param_levels_arg = category_params,
          explanation_text = paste(sm$explanation, "(Excluding reps with divergences.)")
        )
        if (!is.null(plt)) print(plt)
      }
    }
  }
}

# Alpha scatter: avg_est vs. avg_rel_bias
cat("  -> Alpha-parameter scatter plot...\n")
alpha_agg_df <- results_joined_df %>%
  filter(parameter %in% c("alpha[1]", "alpha[2]")) %>%
  group_by(
    T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
    fit_type, fitted_model_code, parameter
  ) %>%
  summarize(
    avg_est      = mean(post_mean, na.rm = TRUE),
    avg_rel_bias = mean(rel_bias, na.rm = TRUE),
    .groups      = "drop"
  )

if (nrow(alpha_agg_df) > 0) {
  alpha_plot <- ggplot(alpha_agg_df, aes(x = avg_est, y = avg_rel_bias, color = fit_type)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    scale_color_viridis_d(option = "plasma", end = 0.8, name = "Fit Type") +
    labs(
      title      = "Alpha Parameter Estimation",
      subtitle   = "Scatter of avg_est (x) vs. avg_rel_bias (y) across conditions",
      x          = "Average Estimated Value (across replications)",
      y          = "Average Relative Bias"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position   = "bottom",
      plot.title        = element_text(hjust = 0.5, size = 12),
      strip.background  = element_blank(),
      panel.grid.major  = element_line(color = "grey90")
    )
  print(alpha_plot)
}

dev.off()
cat("Saved special PDF:", pdf_special, "\n")
cat("\n--- Parameter Recovery Analysis Script Finished ---\n")
