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
DATA_DIR <- file.path(RESULTS_DIR, "../data") # Simplified path
PLOTS_DIR <- file.path(RESULTS_DIR, "plots_param_recovery") # Simplified path
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
cat("Created directory for plots:", PLOTS_DIR, "\n")
data_dir_abs <- normalizePath(DATA_DIR, mustWork = FALSE)
if (!dir.exists(data_dir_abs)) stop("Data directory not found: ", data_dir_abs)

# --- Options ---
`%||%` <- function(a, b) if (!is.null(a)) a else b
theme_set(theme_bw(base_size = 10))
cat("Setup complete.\n")

# --- Load Data ---
sim_conds_file <- file.path(data_dir_abs, "sim_conditions.rds") # Simplified name
if (!file.exists(sim_conds_file)) stop("Sim conditions file not found: ", sim_conds_file)
sim_conditions_df <- readRDS(sim_conds_file) %>%
  mutate(condition_id = as.integer(condition_id)) %>%
  select(condition_id, dgp_copula_type, dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22, correct_fit_model_code) %>%
  distinct(condition_id, .keep_all = TRUE)

param_results_file <- file.path(RESULTS_DIR, "parameter_summary.rds") # Simplified name
if (!file.exists(param_results_file)) stop("Parameter results file not found: ", param_results_file)
results_df <- readRDS(param_results_file) %>% mutate(condition_id = as.integer(condition_id))

sampler_info_file <- file.path(RESULTS_DIR, "sampler_summary.rds") # Simplified name
# Corrected if/else structure for loading sampler info
if (!file.exists(sampler_info_file)) {
  warning("Sampler info file not found.")
  sampler_info_df <- data.frame()
} else {
  sampler_info_df <- readRDS(sampler_info_file) %>% mutate(condition_id = as.integer(condition_id))
} # Corrected structure

if (nrow(results_df) == 0) stop("Parameter results data frame is empty.")
cat("Data loaded successfully.\n")

# --- Define Factor Levels and Order ---
master_param_order <- c("phi11", "phi12", "phi21", "phi22", "mu[1]", "mu[2]", "rho", "theta", "tau", "sigma[1]", "sigma[2]", "xi[1]", "xi[2]", "omega[1]", "omega[2]", "alpha[1]", "alpha[2]")
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
      toupper(substr(ifelse(abs(dgp_alpha1) < 1e-6 & abs(dgp_alpha2) < 1e-6, "normal", "skewnormal"), 1, 1)),
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
    dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)), levels = paste0("Cop=", str_to_title(dgp_copula_levels))),
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
      dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)), levels = paste0("Cop=", str_to_title(dgp_copula_levels))),
      dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
    )
}
cat("Data preparation complete.\n")

# --- Helper Plotting Functions ---
plot_boxplot_dist_agg <- function(data, y_var, category_name, title_suffix, y_label, hline_val = NULL, param_levels_arg) {
  y_sym <- sym(y_var)
  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))
  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
    return(NULL)
  }
  p <- ggplot(data_plot, aes(x = parameter, y = !!y_sym, fill = fit_type)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = ".", outlier.alpha = 0.5, alpha = 0.7, width = 0.7, na.rm = TRUE) +
    {
      if (!is.null(hline_val)) geom_hline(yintercept = hline_val, color = "red", linetype = "dashed")
    } +
    scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fit Type") +
    scale_x_discrete(limits = param_levels_arg, drop = FALSE) +
    labs(title = paste(category_name, "-", title_suffix), x = NULL, y = y_label) +
    # Facet by T/Tau vs Alphas + Copula
    facet_grid(T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac, scales = "free_y", labeller = label_value) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7), legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10), strip.background = element_blank(), panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.3, "lines"), strip.text = element_text(size = 7)
    )
  return(p)
}

plot_rmse_bar_agg <- function(data, category_name, title_suffix, param_levels_arg) {
  data_plot <- data %>%
    filter(param_category == category_name) %>%
    mutate(parameter = factor(parameter, levels = param_levels_arg)) %>%
    filter(!is.na(parameter))
  if (nrow(data_plot) == 0 || length(param_levels_arg) == 0) {
    return(NULL)
  }
  p <- ggplot(data_plot, aes(x = parameter, y = rmse, fill = fit_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30", alpha = 0.8, na.rm = TRUE) +
    scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Fit Type") +
    scale_x_discrete(limits = param_levels_arg, drop = FALSE) +
    labs(title = paste(category_name, "-", title_suffix), x = NULL, y = "RMSE") +
    # Facet by T/Tau vs Alphas + Copula
    facet_grid(T_fac + dgp_tau_fac ~ dgp_alpha1_fac + dgp_alpha2_fac + dgp_copula_fac, scales = "free_y", labeller = label_value) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7), legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10), strip.background = element_blank(), panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.3, "lines"), strip.text = element_text(size = 7)
    )
  return(p)
}
cat("Helper functions defined.\n")

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
    n_reps_param = n(), true_value = first(true_value), avg_est = mean(post_mean, na.rm = TRUE), sd_est = sd(post_mean, na.rm = TRUE),
    avg_bias = mean(bias, na.rm = TRUE), med_bias = median(bias, na.rm = TRUE), avg_rel_bias = mean(rel_bias, na.rm = TRUE),
    med_rel_bias = median(rel_bias, na.rm = TRUE), rmse = sqrt(mean(bias^2, na.rm = TRUE)), mae = mean(abs(bias), na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE), avg_post_sd = mean(post_sd, na.rm = TRUE), avg_ci_width = mean(ci_width, na.rm = TRUE),
    avg_Rhat = mean(Rhat, na.rm = TRUE), pct_rhat_high = mean(Rhat > 1.05, na.rm = TRUE) * 100, med_n_eff = median(n_eff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(parameter = factor(parameter, levels = param_levels), param_category = factor(param_category, levels = category_order))

# --- Join Sampler Summary ---
if (nrow(sampler_info_df) > 0) {
  sampler_summary_agg <- sampler_joined_df %>%
    filter(!is.na(fitted_model_code) & !is.na(dgp_alpha1_fac) & !is.na(dgp_alpha2_fac) & !is.na(dgp_tau_fac)) %>%
    group_by(T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac, fit_type, fitted_model_code) %>%
    summarize(
      n_reps_sampler = n(), avg_divergences = mean(divergences, na.rm = TRUE), pct_reps_diverged = mean(divergences > 0, na.rm = TRUE) * 100,
      avg_maxdepth_exceeded = mean(maxdepth_exceeded, na.rm = TRUE), avg_eFMI = mean(eFMI, na.rm = TRUE), .groups = "drop"
    )
  summary_df <- left_join(summary_df, sampler_summary_agg, by = c("T_fac", "dgp_alpha1_fac", "dgp_alpha2_fac", "dgp_copula_fac", "dgp_tau_fac", "fit_type", "fitted_model_code"))
} else {
  summary_df <- summary_df %>% mutate(n_reps_sampler = NA_integer_, avg_divergences = NA_real_, pct_reps_diverged = NA_real_, avg_maxdepth_exceeded = NA_real_, avg_eFMI = NA_real_)
}

# --- Save Summary Table ---
summary_output_rds <- file.path(RESULTS_DIR, "summary_aggregated.rds") # Simplified name
summary_output_csv <- file.path(RESULTS_DIR, "summary_aggregated.csv") # Simplified name
saveRDS(summary_df, summary_output_rds)
write.csv(summary_df, summary_output_csv, row.names = FALSE)
cat("Saved aggregated summary statistics.\n")

# --- Generate and Save Plots ---
cat("\nGenerating and saving plots...\n")

# --- Aggregated Facet Plots Across Conditions ---
condition_grid_agg <- results_joined_df %>%
  filter(!is.na(T_fac) & !is.na(dgp_alpha1_fac) & !is.na(dgp_alpha2_fac) & !is.na(dgp_copula_fac) & !is.na(dgp_tau_fac)) %>%
  distinct(dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac) %>%
  arrange(dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac)

if (nrow(condition_grid_agg) == 0) {
  warning("No valid conditions found for plotting.", call. = FALSE)
} else {
  cat(paste("Generating plots for", nrow(condition_grid_agg), "DGP Alpha/Copula combinations...\n"))
  plot_count <- 0
  purrr::pwalk(condition_grid_agg, function(dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac) {
    a1_label <- gsub("[^a-zA-Z0-9-]", "_", dgp_alpha1_fac)
    a2_label <- gsub("[^a-zA-Z0-9-]", "_", dgp_alpha2_fac)
    cop_label <- gsub("[^a-zA-Z0-9]", "_", dgp_copula_fac)
    cond_title_base <- sprintf("DGP: %s/%s, %s", dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac)

    current_results_subset <- results_joined_df %>%
      filter(dgp_alpha1_fac == !!dgp_alpha1_fac,
             dgp_alpha2_fac == !!dgp_alpha2_fac,
             dgp_copula_fac == !!dgp_copula_fac,
             !is.na(dgp_tau_fac))
    current_summary_subset <- summary_df %>%
      filter(dgp_alpha1_fac == !!dgp_alpha1_fac,
             dgp_alpha2_fac == !!dgp_alpha2_fac,
             dgp_copula_fac == !!dgp_copula_fac,
             !is.na(dgp_tau_fac))
    if (nrow(current_results_subset) == 0) return()

    for (cat_name in category_order) {
      cat_label <- gsub("[^a-zA-Z0-9]", "_", gsub("/", "_", cat_name))
      category_params <- intersect(master_param_order,
                                   unique(as.character(current_results_subset %>% filter(param_category == cat_name) %>% pull(parameter))))
      if (length(category_params) == 0 || nrow(current_results_subset %>% filter(param_category == cat_name)) == 0) next

      plot_types <- list(
        list(metric = "bias", label = "Bias", hline = 0),
        list(metric = "rel_bias_capped", label = "RelBias", hline = 0),
        list(metric = "post_sd", label = "PostSD", hline = NULL),
        list(metric = "rmse", label = "RMSE", hline = NULL, type = "bar")
      )
      for (pt in plot_types) {
        plot_filename <- file.path(PLOTS_DIR, sprintf("plot_%s_%s_%s_%s_%s.png",
                                                   pt$label, cat_label, a1_label, a2_label, cop_label))
        plot_obj <- NULL
        tryCatch({
          if (pt$type %||% "box" == "box") {
            plot_data <- current_results_subset
            if (pt$metric == "rel_bias_capped") {
              plot_data <- plot_data %>%
                filter(abs(true_value) > 1e-6 & !is.na(true_value) & is.finite(bias)) %>%
                mutate(rel_bias_capped = pmax(-2, pmin(2, rel_bias)))
              if (nrow(plot_data) == 0) stop("No data for rel bias")
            }
            plot_obj <- plot_boxplot_dist_agg(plot_data, pt$metric, cat_name, cond_title_base,
                                              pt$label, pt$hline, category_params)
          } else if (pt$type == "bar") {
            if (nrow(current_summary_subset %>% filter(param_category == cat_name)) > 0) {
              plot_obj <- plot_rmse_bar_agg(current_summary_subset, cat_name, cond_title_base, category_params)
            } else {
              stop("No summary data for RMSE")
            }
          }
          if (!is.null(plot_obj)) {
            ggsave(filename = plot_filename, plot = plot_obj, width = 11, height = 6, units = "in", dpi = 150)
            plot_count <- plot_count + 1
          }
        }, error = function(e) {
          cat(paste0("     ERROR plot ", basename(plot_filename), ": ", e$message, "\n"))
        })
      }
    }
  })
  cat("\nFinished plots. Total saved:", plot_count, ". Plots in:", PLOTS_DIR, "\n")
}

# --- Per-Condition PDF Plots ---
condition_info_df <- sim_conditions_df %>%
  mutate(
    T_fac = factor(paste0("T=", T), levels = paste0("T=", T_levels)),
    dgp_alpha1_fac = factor(paste0("A1=", dgp_alpha1), levels = paste0("A1=", dgp_alpha1_levels)),
    dgp_alpha2_fac = factor(paste0("A2=", dgp_alpha2), levels = paste0("A2=", dgp_alpha2_levels)),
    dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)), levels = paste0("Cop=", str_to_title(dgp_copula_levels))),
    dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
  ) %>%
  select(condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
         phi11, phi12, phi21, phi22)

cat("Generating per-condition parameter recovery PDFs...\n")
purrr::pwalk(condition_info_df,
             function(condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
                      phi11, phi12, phi21, phi22) {
  cond_res <- results_joined_df %>% filter(condition_id == !!condition_id)
  cond_sum <- summary_df %>% filter(condition_id == !!condition_id)
  if (nrow(cond_res) == 0) return()
  pdf_file <- file.path(PLOTS_DIR, sprintf("param_recov_cond_%03d.pdf", condition_id))
  grDevices::pdf(pdf_file, width = 11, height = 6)
  grid::grid.newpage()
  grid::grid.text(sprintf("Condition %03d\n%s %s %s %s %s\nphi11=%.2f phi12=%.2f phi21=%.2f phi22=%.2f",
                          condition_id, T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac,
                          phi11, phi12, phi21, phi22), gp = grid::gpar(cex = 0.9))
  plot_types <- list(
    list(metric = "bias", label = "Bias", hline = 0),
    list(metric = "rel_bias_capped", label = "RelBias", hline = 0),
    list(metric = "post_sd", label = "PostSD", hline = NULL),
    list(metric = "rmse", label = "RMSE", hline = NULL, type = "bar")
  )
  for (cat_name in category_order) {
    category_params <- intersect(master_param_order,
                                 unique(as.character(cond_res %>% filter(param_category == cat_name) %>% pull(parameter))))
    if (length(category_params) == 0 || nrow(cond_res %>% filter(param_category == cat_name)) == 0) next
    for (pt in plot_types) {
      plot_obj <- NULL
      tryCatch({
        if (pt$type %||% "box" == "box") {
          plot_data <- cond_res
          if (pt$metric == "rel_bias_capped") {
            plot_data <- plot_data %>%
              filter(abs(true_value) > 1e-6 & !is.na(true_value) & is.finite(bias)) %>%
              mutate(rel_bias_capped = pmax(-2, pmin(2, rel_bias)))
            if (nrow(plot_data) == 0) stop("No data for rel bias")
          }
          plot_obj <- plot_boxplot_dist_agg(plot_data, pt$metric, cat_name,
                                            paste(T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac),
                                            pt$label, pt$hline, category_params)
        } else if (pt$type == "bar") {
          if (nrow(cond_sum %>% filter(param_category == cat_name)) > 0) {
            plot_obj <- plot_rmse_bar_agg(cond_sum, cat_name,
                                         paste(T_fac, dgp_alpha1_fac, dgp_alpha2_fac, dgp_copula_fac, dgp_tau_fac),
                                         category_params)
          } else {
            stop("No summary data for RMSE")
          }
        }
        if (!is.null(plot_obj)) print(plot_obj)
      }, error = function(e) {
        cat(paste0("     ERROR plot cond", condition_id, ": ", e$message, "\n"))
      })
    }
  }
  grDevices::dev.off()
cat("Saved", pdf_file, "\n")
})
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\n--- Parameter Recovery Analysis Script Finished ---\n")
