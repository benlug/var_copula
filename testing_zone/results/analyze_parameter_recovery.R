###########################################################################
# analyze_parameter_recovery.R
#
# Analyzes parameter recovery from VAR(1) Copula simulation (separate alphas).
# Simplified naming and plotting structure.
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
RESULTS_DIR <- getwd()
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
}
cat("Data preparation complete.\n")

# --- Helper Plotting Functions ---
plot_boxplot_dist_agg <- function(data, y_var, category_name, title_suffix, y_label, hline_val = NULL, param_levels_arg, explanation_text = "") {
  # This function creates a boxplot over replications for metric y_var,
  # faceted by T/tau vs alpha/copula, with param_category filtering.
  # 'explanation_text' is added to the plot subtitle for interpretation hints.

  y_sym <- sym(y_var)
  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))

  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
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
    # Facet by T/Tau vs Alphas + Copula
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
  # This function creates a bar chart of RMSE (aggregated across replications).

  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))

  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
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
      title = paste(category_name, "-", title_suffix),
      subtitle = explanation_text,
      x = NULL,
      y = "RMSE"
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

# --- Calculate Summary Statistics ---
cat("Calculating aggregated summary statistics...\n")
summary_df <- results_joined_df %>%
  filter(!is.na(parameter) & !is.na(fitted_model_code) & !is.na(param_category)) %>%
  group_by(
    T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
    fit_type, fitted_model_code, parameter, param_category
  ) %>%
  filter(!is.na(dgp_alpha1_fac) & !is.na(dgp_alpha2_fac) & !is.na(dgp_tau_fac)) %>%
  summarize(
    n_reps_param   = n(),
    true_value     = first(true_value),
    avg_est        = mean(post_mean, na.rm = TRUE),
    sd_est         = sd(post_mean, na.rm = TRUE),
    avg_bias       = mean(bias, na.rm = TRUE),
    med_bias       = median(bias, na.rm = TRUE),
    avg_rel_bias   = mean(rel_bias, na.rm = TRUE),
    med_rel_bias   = median(rel_bias, na.rm = TRUE),
    rmse           = sqrt(mean(bias^2, na.rm = TRUE)),
    mae            = mean(abs(bias), na.rm = TRUE),
    coverage       = mean(coverage, na.rm = TRUE),
    avg_post_sd    = mean(post_sd, na.rm = TRUE),
    avg_ci_width   = mean(ci_width, na.rm = TRUE),
    avg_Rhat       = mean(Rhat, na.rm = TRUE),
    pct_rhat_high  = mean(Rhat > 1.05, na.rm = TRUE) * 100,
    med_n_eff      = median(n_eff, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(
    parameter      = factor(parameter, levels = param_levels),
    param_category = factor(param_category, levels = category_order)
  )

# --- Join Sampler Summary (optional) ---
if (nrow(sampler_info_df) > 0) {
  sampler_summary_agg <- sampler_joined_df %>%
    filter(!is.na(fitted_model_code) & !is.na(dgp_alpha1_fac) & !is.na(dgp_alpha2_fac) & !is.na(dgp_tau_fac)) %>%
    group_by(T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, fit_type, fitted_model_code) %>%
    summarize(
      n_reps_sampler      = n(),
      avg_divergences     = mean(divergences, na.rm = TRUE),
      pct_reps_diverged   = mean(divergences > 0, na.rm = TRUE) * 100,
      avg_maxdepth_exceed = mean(maxdepth_exceeded, na.rm = TRUE),
      avg_eFMI            = mean(eFMI, na.rm = TRUE),
      .groups             = "drop"
    )
  summary_df <- left_join(
    summary_df,
    sampler_summary_agg,
    by = c("T_fac", "dgp_alpha1_fac", "dgp_alpha2_fac", "dgp_copula_fac", "dgp_tau_fac", "fit_type", "fitted_model_code")
  )
} else {
  summary_df <- summary_df %>%
    mutate(
      n_reps_sampler      = NA_integer_,
      avg_divergences     = NA_real_,
      pct_reps_diverged   = NA_real_,
      avg_maxdepth_exceed = NA_real_,
      avg_eFMI            = NA_real_
    )
}

# --- Save Summary Table ---
summary_output_rds <- file.path(RESULTS_DIR, "summary_aggregated.rds")
summary_output_csv <- file.path(RESULTS_DIR, "summary_aggregated.csv")
saveRDS(summary_df, summary_output_rds)
write.csv(summary_df, summary_output_csv, row.names = FALSE)
cat("Saved aggregated summary statistics.\n")

# -----------------------------------------------------------------------------
# Generate Plots
# We now separate them into two PDF files: main vs. detailed
# -----------------------------------------------------------------------------
cat("\nGenerating faceted parameter recovery plots...\n")

# 1) MAIN PDF: relative bias + postSD
pdf_main <- file.path(PLOTS_DIR, "parameter_recovery_faceted_main.pdf")
grDevices::pdf(pdf_main, width = 11, height = 6)
cat("Creating MAIN PDF with relative bias + posterior SD...\n")

# We'll define the two main metrics with short explanations:
main_plots <- list(
  list(
    metric         = "rel_bias_capped",
    label          = "Relative Bias",
    hline          = 0,
    explanation    = "Relative bias = (Estimate - True) / |True|. 0 means no bias. Positive = overestimate, negative = underestimate.",
    type           = "box"
  ),
  list(
    metric         = "post_sd",
    label          = "Posterior SD",
    explanation    = "Posterior SD = average standard deviation of the posterior for each parameter. Higher = more uncertainty.",
    type           = "box"
  )
)

# We'll create 'rel_bias_capped' on-the-fly in results_joined_df
# limiting rel_bias to [-2, 2] for visibility:
results_joined_df_main <- results_joined_df %>%
  mutate(
    rel_bias_capped = ifelse(
      !is.na(rel_bias),
      pmax(-2, pmin(2, rel_bias)),
      NA_real_
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
    tryCatch(
      {
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
      },
      error = function(e) {
        cat("ERROR in main metric plotting: ", e$message, "\n")
      }
    )
  }
}

dev.off()
cat("Saved main PDF:", pdf_main, "\n")

# 2) DETAILED PDF: raw bias + RMSE
pdf_detail <- file.path(PLOTS_DIR, "parameter_recovery_faceted_detailed.pdf")
grDevices::pdf(pdf_detail, width = 11, height = 6)
cat("Creating DETAILED PDF with raw bias + RMSE...\n")

detailed_plots <- list(
  list(
    metric         = "bias",
    label          = "Raw Bias",
    hline          = 0,
    explanation    = "Raw bias = Estimate - True. 0 means no bias. Positive = overestimate, negative = underestimate.",
    type           = "box"
  ),
  list(
    metric         = "rmse",
    label          = "RMSE",
    explanation    = "Root Mean Squared Error across replications. Lower = better recovery.",
    type           = "bar"
  )
)

for (cat_name in category_order) {
  category_params <- intersect(
    master_param_order,
    unique(as.character(results_joined_df %>% filter(param_category == cat_name) %>% pull(parameter)))
  )
  if (length(category_params) == 0) next

  for (plot_info in detailed_plots) {
    plot_obj <- NULL
    tryCatch(
      {
        if (plot_info$type == "box") {
          # For raw bias
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
        } else if (plot_info$type == "bar") {
          # For RMSE (aggregated in summary_df)
          if (nrow(summary_df %>% filter(param_category == cat_name)) > 0) {
            plot_obj <- plot_rmse_bar_agg(
              data = summary_df %>% filter(param_category == cat_name),
              category_name = cat_name,
              title_suffix = "All Conditions",
              param_levels_arg = category_params,
              explanation_text = plot_info$explanation
            )
          }
        }
        if (!is.null(plot_obj)) print(plot_obj)
      },
      error = function(e) {
        cat("ERROR in detailed metric plotting: ", e$message, "\n")
      }
    )
  }
}

dev.off()
cat("Saved detailed PDF:", pdf_detail, "\n")

cat("\n--- Parameter Recovery Analysis Script Finished ---\n")
