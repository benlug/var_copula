###########################################################################
# analyze_convergence.R
#
# Comprehensive MCMC convergence analysis for VAR(1) Copula simulation
# Creates detailed PDF reports for each condition with tables and plots
#
# This script helps us understand whether our Bayesian inference via Stan
# is working properly. Think of it as a health check for our statistical
# machinery - we're making sure the engine is running smoothly before
# we trust the results it produces.
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
library(gridExtra)

# --- Configuration ---
RESULTS_DIR <- getwd()
DATA_DIR <- file.path(RESULTS_DIR, "../data")
PLOTS_DIR <- file.path(RESULTS_DIR, "convergence_reports")
FITS_DIR <- file.path(RESULTS_DIR, "../fits")

# Create output directory for our diagnostic reports
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
cat("Created directory for convergence reports:", PLOTS_DIR, "\n")

# --- Options ---
`%||%` <- function(a, b) if (!is.null(a)) a else b
theme_set(theme_bw(base_size = 10))

# Set up consistent color scheme for our diagnostic plots
diagnostic_colors <- list(
  good = "#2E7D32", # Green for good values
  warning = "#F57C00", # Orange for warning values
  bad = "#D32F2F", # Red for problematic values
  threshold = "#1976D2" # Blue for threshold lines
)

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

#' Read Stan fit safely with error handling
#'
#' This function protects us from corrupted files or reading errors
#' that could crash our entire analysis pipeline.
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

#' Calculate Energy Fraction of Missing Information (E-FMI)
#'
#' E-FMI tells us how well the momentum resampling in HMC is working.
#' Think of it as measuring how efficiently our sampler is exploring
#' the parameter space. Low values suggest the sampler is struggling.
calc_e_fmi <- function(energy_vec) {
  n <- length(energy_vec)
  if (n < 3 || any(!is.finite(energy_vec))) {
    return(NA_real_)
  }

  # Calculate variance of energy
  E_mean <- mean(energy_vec)
  var_E <- var(energy_vec)

  # Calculate lag-1 autocorrelation of energy
  lag1_sum <- sum((energy_vec[-1] - E_mean) * (energy_vec[-n] - E_mean))
  lag1_cov <- lag1_sum / (n - 1)

  # E-FMI formula: Var(E) / (Var(E) + 2*Cov(E_t, E_{t-1}))
  denom <- var_E + 2 * lag1_cov
  if (is.na(denom) || denom <= 1e-9) {
    return(NA_real_)
  }

  efmi_val <- var_E / denom
  return(pmax(0, pmin(1.5, efmi_val))) # Clamp to reasonable range
}

#' Create diagnostic status indicator
#'
#' This helper function creates a colored text indicator showing
#' whether a diagnostic value is good, warning, or bad.
diagnostic_status <- function(value, type = "divergences", is_percentage = FALSE) {
  if (is.na(value)) {
    return("N/A")
  }

  # Define thresholds for each diagnostic type
  thresholds <- switch(type,
    "divergences" = list(good = 0, warning = 1, bad = 10),
    "efmi" = list(good = 0.3, warning = 0.2, bad = 0.1),
    "rhat" = list(good = 1.01, warning = 1.05, bad = 1.1),
    "neff_ratio" = list(good = 0.1, warning = 0.05, bad = 0.01),
    "accept_stat" = list(good = 0.9, warning = 0.8, bad = 0.7),
    list(good = NA, warning = NA, bad = NA)
  )

  # Format the value
  formatted_value <- if (is_percentage) {
    sprintf("%.1f%%", value)
  } else if (type == "divergences") {
    sprintf("%.0f", value)
  } else {
    sprintf("%.3f", value)
  }

  # Determine status and color
  if (type %in% c("divergences", "rhat")) {
    # For these, lower is better
    if (value <= thresholds$good) {
      status <- "GOOD"
      color <- diagnostic_colors$good
    } else if (value <= thresholds$warning) {
      status <- "WARNING"
      color <- diagnostic_colors$warning
    } else {
      status <- "BAD"
      color <- diagnostic_colors$bad
    }
  } else {
    # For these, higher is better
    if (value >= thresholds$good) {
      status <- "GOOD"
      color <- diagnostic_colors$good
    } else if (value >= thresholds$warning) {
      status <- "WARNING"
      color <- diagnostic_colors$warning
    } else {
      status <- "BAD"
      color <- diagnostic_colors$bad
    }
  }

  return(sprintf("%s [%s]", formatted_value, status))
}

#' Create a summary table for a single condition
#'
#' This function creates a nicely formatted table summarizing all
#' convergence diagnostics for a given condition.
create_condition_summary_table <- function(condition_data, sampler_data) {
  # Calculate key metrics for each model type fitted
  summary_by_model <- sampler_data %>%
    group_by(fitted_model_code, fit_type) %>%
    summarise(
      n_fits = n(),
      # Divergences
      n_diverged = sum(divergences > 0, na.rm = TRUE),
      pct_diverged = mean(divergences > 0, na.rm = TRUE) * 100,
      avg_divergences = mean(divergences[divergences > 0], na.rm = TRUE),
      max_divergences = max(divergences, na.rm = TRUE),
      # E-FMI
      avg_efmi = mean(eFMI, na.rm = TRUE),
      min_efmi = min(eFMI, na.rm = TRUE),
      n_low_efmi = sum(eFMI < 0.3, na.rm = TRUE),
      # Tree depth
      n_max_treedepth = sum(maxdepth_exceeded > 0, na.rm = TRUE),
      # Acceptance statistic
      avg_accept_stat = mean(avg_accept_stat, na.rm = TRUE),
      .groups = "drop"
    )

  # Add R-hat and N_eff information if available
  if (nrow(condition_data) > 0) {
    param_summary <- condition_data %>%
      group_by(fitted_model_code, fit_type) %>%
      summarise(
        avg_rhat = mean(Rhat, na.rm = TRUE),
        max_rhat = max(Rhat, na.rm = TRUE),
        n_high_rhat = sum(Rhat > 1.05, na.rm = TRUE),
        median_neff = median(n_eff, na.rm = TRUE),
        min_neff = min(n_eff, na.rm = TRUE),
        .groups = "drop"
      )

    summary_by_model <- left_join(summary_by_model, param_summary,
      by = c("fitted_model_code", "fit_type")
    )
  }

  # Create formatted table
  summary_table <- summary_by_model %>%
    mutate(
      Model = paste0(fitted_model_code, " (", fit_type, ")"),
      Divergences = diagnostic_status(avg_divergences, "divergences"),
      `Max Divergences` = diagnostic_status(max_divergences, "divergences"),
      `% with Divergences` = diagnostic_status(pct_diverged, "divergences", TRUE),
      `E-FMI` = diagnostic_status(avg_efmi, "efmi"),
      `Min E-FMI` = diagnostic_status(min_efmi, "efmi"),
      `Avg R-hat` = diagnostic_status(avg_rhat, "rhat"),
      `Max R-hat` = diagnostic_status(max_rhat, "rhat"),
      `Median N_eff` = sprintf("%.0f", median_neff),
      `Accept Stat` = diagnostic_status(avg_accept_stat, "accept_stat")
    ) %>%
    select(
      Model, Divergences, `% with Divergences`, `E-FMI`, `Min E-FMI`,
      `Avg R-hat`, `Max R-hat`, `Median N_eff`, `Accept Stat`
    )

  return(summary_table)
}

#' Create parameter-specific diagnostic table
#'
#' This breaks down convergence by parameter type, helping us identify
#' which parameters are causing sampling difficulties.
create_parameter_diagnostic_table <- function(condition_data) {
  if (nrow(condition_data) == 0) {
    return(data.frame(Parameter = "No data available"))
  }

  param_summary <- condition_data %>%
    group_by(param_category, fitted_model_code, fit_type) %>%
    summarise(
      n_params = n_distinct(parameter),
      avg_rhat = mean(Rhat, na.rm = TRUE),
      pct_high_rhat = mean(Rhat > 1.05, na.rm = TRUE) * 100,
      median_neff = median(n_eff, na.rm = TRUE),
      min_neff = min(n_eff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      `Parameter Group` = param_category,
      Model = paste0(fitted_model_code, " (", fit_type, ")"),
      `Avg R-hat` = diagnostic_status(avg_rhat, "rhat"),
      `% R-hat > 1.05` = diagnostic_status(pct_high_rhat, "divergences", TRUE),
      `Median N_eff` = sprintf("%.0f", median_neff),
      `Min N_eff` = sprintf("%.0f", min_neff)
    ) %>%
    select(
      `Parameter Group`, Model, `Avg R-hat`, `% R-hat > 1.05`,
      `Median N_eff`, `Min N_eff`
    ) %>%
    arrange(`Parameter Group`, Model)

  return(param_summary)
}

#' Create diagnostic plots for a condition
#'
#' This function generates a comprehensive set of plots to visualize
#' convergence diagnostics. Each plot tells us something different
#' about how well our sampler is working.
create_diagnostic_plots <- function(condition_data, sampler_data, condition_info) {
  plots <- list()

  # Common plot styling
  base_theme <- theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 9)
    )

  # 1. Divergences plot
  if (nrow(sampler_data) > 0) {
    p_div <- ggplot(sampler_data, aes(x = fit_type, y = divergences, fill = fitted_model_code)) +
      geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
      geom_boxplot(
        width = 0.15, position = position_dodge(width = 0.9),
        outlier.shape = 21, outlier.alpha = 0.5, alpha = 0.8
      ) +
      scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Model") +
      scale_y_sqrt(breaks = c(0, 1, 5, 10, 25, 50, 100)) + # Square root scale for better visibility
      labs(
        title = "Divergent Transitions by Model Type",
        subtitle = "Lower is better - ideally zero divergences",
        x = "Model Specification",
        y = "Number of Divergences (sqrt scale)"
      ) +
      base_theme

    plots$divergences <- p_div

    # 2. E-FMI plot
    p_efmi <- ggplot(sampler_data, aes(x = fit_type, y = eFMI, fill = fitted_model_code)) +
      geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
      geom_boxplot(
        width = 0.15, position = position_dodge(width = 0.9),
        outlier.shape = 21, outlier.alpha = 0.5, alpha = 0.8
      ) +
      geom_hline(
        yintercept = 0.3, linetype = "dashed",
        color = diagnostic_colors$threshold, size = 1
      ) +
      scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Model") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Energy Fraction of Missing Information (E-FMI)",
        subtitle = "Should be above 0.3 (blue line) for good sampling",
        x = "Model Specification",
        y = "E-FMI"
      ) +
      base_theme

    plots$efmi <- p_efmi

    # 3. Acceptance statistic plot
    p_accept <- ggplot(sampler_data, aes(x = fit_type, y = avg_accept_stat, fill = fitted_model_code)) +
      geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
      geom_boxplot(
        width = 0.15, position = position_dodge(width = 0.9),
        outlier.shape = 21, outlier.alpha = 0.5, alpha = 0.8
      ) +
      geom_hline(
        yintercept = 0.99, linetype = "dashed",
        color = diagnostic_colors$threshold, size = 1
      ) +
      scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Model") +
      scale_y_continuous(limits = c(0.7, 1)) +
      labs(
        title = "Average Acceptance Statistic",
        subtitle = "Target is 0.99 (blue line) based on adapt_delta setting",
        x = "Model Specification",
        y = "Acceptance Probability"
      ) +
      base_theme

    plots$accept_stat <- p_accept
  }

  # 4. R-hat distribution plot
  if (nrow(condition_data) > 0) {
    # Calculate average R-hat per fit
    avg_rhat_data <- condition_data %>%
      group_by(rep_i, fitted_model_code, fit_type) %>%
      summarise(avg_Rhat = mean(Rhat, na.rm = TRUE), .groups = "drop")

    p_rhat <- ggplot(avg_rhat_data, aes(x = fit_type, y = avg_Rhat, fill = fitted_model_code)) +
      geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
      geom_boxplot(
        width = 0.15, position = position_dodge(width = 0.9),
        outlier.shape = 21, outlier.alpha = 0.5, alpha = 0.8
      ) +
      geom_hline(
        yintercept = 1.05, linetype = "dashed",
        color = diagnostic_colors$threshold, size = 1
      ) +
      scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Model") +
      coord_cartesian(ylim = c(0.99, 1.1)) +
      labs(
        title = "Potential Scale Reduction Factor (R-hat)",
        subtitle = "Should be below 1.05 (blue line) for all parameters",
        x = "Model Specification",
        y = "Average R-hat"
      ) +
      base_theme

    plots$rhat <- p_rhat

    # 5. Effective sample size plot
    med_neff_data <- condition_data %>%
      group_by(rep_i, fitted_model_code, fit_type) %>%
      summarise(med_n_eff = median(n_eff, na.rm = TRUE), .groups = "drop")

    p_neff <- ggplot(med_neff_data, aes(x = fit_type, y = med_n_eff, fill = fitted_model_code)) +
      geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
      geom_boxplot(
        width = 0.15, position = position_dodge(width = 0.9),
        outlier.shape = 21, outlier.alpha = 0.5, alpha = 0.8
      ) +
      geom_hline(
        yintercept = 400, linetype = "dashed",
        color = diagnostic_colors$threshold, size = 1
      ) +
      scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Model") +
      scale_y_log10(limits = c(10, NA), breaks = c(10, 100, 400, 1000, 4000)) +
      annotation_logticks(sides = "l") +
      labs(
        title = "Effective Sample Size",
        subtitle = "Higher is better - ideally above 400 (blue line) for key parameters",
        x = "Model Specification",
        y = "Median N_eff (log scale)"
      ) +
      base_theme

    plots$neff <- p_neff

    # 6. Parameter-specific R-hat plot
    p_param_rhat <- ggplot(condition_data, aes(x = parameter, y = Rhat, color = fit_type)) +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      geom_hline(
        yintercept = 1.05, linetype = "dashed",
        color = diagnostic_colors$threshold, size = 1
      ) +
      facet_wrap(~param_category, scales = "free_x", ncol = 2) +
      scale_color_viridis_d(option = "plasma", end = 0.8, name = "Model Type") +
      coord_cartesian(ylim = c(0.99, 1.1)) +
      labs(
        title = "R-hat by Parameter Type",
        subtitle = "Helps identify which parameters are difficult to sample",
        x = NULL,
        y = "R-hat"
      ) +
      base_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    plots$param_rhat <- p_param_rhat
  }

  return(plots)
}

#' Create example trace plots for problematic fits
#'
#' When we see convergence problems, trace plots help us understand
#' what's going wrong. This function finds the most problematic fit
#' and shows its traces.
create_trace_plot_for_worst_fit <- function(condition_id, sampler_data, fits_dir) {
  # Find the fit with the most divergences
  worst_fit <- sampler_data %>%
    arrange(desc(divergences)) %>%
    slice(1)

  if (nrow(worst_fit) == 0 || worst_fit$divergences == 0) {
    return(NULL)
  }

  # Load the problematic fit
  fit_file <- file.path(
    fits_dir,
    sprintf(
      "fit_%s_cond%03d_rep%03d.rds",
      worst_fit$fitted_model_code,
      condition_id,
      worst_fit$rep_i
    )
  )

  fit <- safe_read_stanfit(fit_file)
  if (is.null(fit)) {
    return(NULL)
  }

  # Extract key parameters for trace plot
  param_names <- c("phi11", "phi12", "phi21", "phi22", "rho", "theta")
  available_params <- param_names[param_names %in% names(fit)]

  if (length(available_params) == 0) {
    return(NULL)
  }

  # Create trace plot using bayesplot
  trace_plot <- tryCatch(
    {
      bayesplot::mcmc_trace(fit, pars = available_params[1:min(4, length(available_params))]) +
        labs(
          title = sprintf(
            "Trace Plot: %s Model (Rep %d with %d divergences)",
            worst_fit$fitted_model_code,
            worst_fit$rep_i,
            worst_fit$divergences
          ),
          subtitle = "Looking for 'hairy caterpillar' appearance - good mixing across chains"
        ) +
        theme_bw(base_size = 10)
    },
    error = function(e) NULL
  )

  return(trace_plot)
}

#' Generate comprehensive PDF report for a condition
#'
#' This is our main reporting function. It creates a complete diagnostic
#' report for a single simulation condition, combining tables and plots
#' into an easy-to-read PDF document.
generate_condition_report <- function(condition_id, condition_info,
                                      condition_data, sampler_data,
                                      output_dir, fits_dir) {
  # Create filename for this condition's report
  pdf_file <- file.path(output_dir, sprintf("convergence_report_cond_%03d.pdf", condition_id))

  # Set up PDF device with appropriate size
  pdf(pdf_file, width = 11, height = 8.5)

  # Page 1: Title and condition information
  grid.newpage()

  # Create title
  title_text <- sprintf("MCMC Convergence Report\nCondition %03d", condition_id)
  condition_text <- sprintf(
    "\nSimulation Parameters:\n%s | %s | %s | %s | %s\nVAR Coefficients: φ11=%.2f, φ12=%.2f, φ21=%.2f, φ22=%.2f",
    condition_info$T_fac,
    condition_info$dgp_alpha1_fac,
    condition_info$dgp_alpha2_fac,
    condition_info$dgp_copula_fac,
    condition_info$dgp_tau_fac,
    condition_info$phi11,
    condition_info$phi12,
    condition_info$phi21,
    condition_info$phi22
  )

  grid.text(title_text,
    x = 0.5, y = 0.8,
    gp = gpar(fontsize = 20, fontface = "bold")
  )
  grid.text(condition_text,
    x = 0.5, y = 0.6,
    gp = gpar(fontsize = 12)
  )

  # Add interpretation guide
  guide_text <- "Interpretation Guide:
• Divergences: Should be 0, indicates numerical problems if present
• E-FMI: Should be > 0.3, measures sampling efficiency
• R-hat: Should be < 1.05, measures chain convergence
• N_eff: Should be > 400 for key parameters, effective sample size
• Accept Stat: Should be near target (0.99), adaptation quality"

  grid.text(guide_text,
    x = 0.5, y = 0.3,
    gp = gpar(fontsize = 10, fontfamily = "mono"),
    just = "center"
  )

  # Page 2: Summary tables
  grid.newpage()
  grid.text("Convergence Summary Tables",
    x = 0.5, y = 0.95,
    gp = gpar(fontsize = 16, fontface = "bold")
  )

  # Create and display overall summary table
  summary_table <- create_condition_summary_table(condition_data, sampler_data)

  if (nrow(summary_table) > 0) {
    # Convert to grob for plotting
    summary_grob <- tableGrob(
      summary_table,
      theme = ttheme_default(
        base_size = 9,
        core = list(fg_params = list(hjust = 0, x = 0.1)),
        colhead = list(bg_params = list(fill = "lightblue"))
      )
    )

    grid.draw(summary_grob)
  }

  # Page 3: Parameter-specific diagnostics table
  if (nrow(condition_data) > 0) {
    grid.newpage()
    grid.text("Parameter-Specific Diagnostics",
      x = 0.5, y = 0.95,
      gp = gpar(fontsize = 16, fontface = "bold")
    )

    param_table <- create_parameter_diagnostic_table(condition_data)

    param_grob <- tableGrob(
      param_table,
      theme = ttheme_default(
        base_size = 8,
        core = list(fg_params = list(hjust = 0, x = 0.1)),
        colhead = list(bg_params = list(fill = "lightgreen"))
      )
    )

    # Position the table
    grid.draw(param_grob)
  }

  # Pages 4+: Diagnostic plots
  plots <- create_diagnostic_plots(condition_data, sampler_data, condition_info)

  # Plot each diagnostic
  for (plot_name in names(plots)) {
    if (!is.null(plots[[plot_name]])) {
      grid.newpage()
      print(plots[[plot_name]])
    }
  }

  # Final page: Trace plot for worst fit (if applicable)
  trace_plot <- create_trace_plot_for_worst_fit(condition_id, sampler_data, fits_dir)
  if (!is.null(trace_plot)) {
    grid.newpage()
    print(trace_plot)
  }

  # Close PDF device
  dev.off()

  cat(sprintf("Generated report for condition %03d: %s\n", condition_id, pdf_file))
}

# ===========================================================================
# MAIN ANALYSIS WORKFLOW
# ===========================================================================

cat("\n=== Starting Comprehensive Convergence Analysis ===\n")

# Step 1: Load simulation conditions
sim_conds_file <- file.path(DATA_DIR, "sim_conditions.rds")
if (!file.exists(sim_conds_file)) {
  stop("Simulation conditions file not found: ", sim_conds_file)
}

sim_conditions_df <- readRDS(sim_conds_file) %>%
  mutate(condition_id = as.integer(condition_id)) %>%
  select(
    condition_id, dgp_copula_type, dgp_alpha1, dgp_alpha2,
    dgp_tau, T, phi11:phi22, correct_fit_model_code
  ) %>%
  distinct(condition_id, .keep_all = TRUE)

cat(sprintf("Loaded %d simulation conditions.\n", nrow(sim_conditions_df)))

# Step 2: Load parameter results
param_results_file <- file.path(RESULTS_DIR, "parameter_summary.rds")
if (!file.exists(param_results_file)) {
  stop("Parameter results file not found: ", param_results_file)
}

results_df <- readRDS(param_results_file) %>%
  mutate(condition_id = as.integer(condition_id))

cat(sprintf("Loaded parameter results: %d records.\n", nrow(results_df)))

# Step 3: Load sampler diagnostics
sampler_info_file <- file.path(RESULTS_DIR, "sampler_summary.rds")
if (!file.exists(sampler_info_file)) {
  warning("Sampler info file not found. Some diagnostics will be unavailable.")
  sampler_info_df <- data.frame()
} else {
  sampler_info_df <- readRDS(sampler_info_file) %>%
    mutate(condition_id = as.integer(condition_id))
  cat(sprintf("Loaded sampler diagnostics: %d records.\n", nrow(sampler_info_df)))
}

# Step 4: Prepare factor levels for consistent plotting
dgp_copula_levels <- sort(unique(c(results_df$dgp_copula_type, sampler_info_df$dgp_copula_type)))
dgp_alpha1_levels <- sort(unique(c(results_df$dgp_alpha1, sampler_info_df$dgp_alpha1)))
dgp_alpha2_levels <- sort(unique(c(results_df$dgp_alpha2, sampler_info_df$dgp_alpha2)))
dgp_tau_levels <- sort(unique(c(results_df$dgp_tau, sampler_info_df$dgp_tau)))
T_levels <- sort(unique(c(results_df$T, sampler_info_df$T)))

# Create condition info with factors
condition_info_df <- sim_conditions_df %>%
  mutate(
    T_fac = factor(paste0("T=", T), levels = paste0("T=", T_levels)),
    dgp_alpha1_fac = factor(paste0("A1=", dgp_alpha1), levels = paste0("A1=", dgp_alpha1_levels)),
    dgp_alpha2_fac = factor(paste0("A2=", dgp_alpha2), levels = paste0("A2=", dgp_alpha2_levels)),
    dgp_copula_fac = factor(paste0("Cop=", str_to_title(dgp_copula_type)),
      levels = paste0("Cop=", str_to_title(dgp_copula_levels))
    ),
    dgp_tau_fac = factor(paste0("Tau=", dgp_tau), levels = paste0("Tau=", dgp_tau_levels))
  )

# Step 5: Generate reports for each condition
cat("\n--- Generating Individual Condition Reports ---\n")

# Process each condition
for (i in 1:nrow(condition_info_df)) {
  condition <- condition_info_df[i, ]
  cond_id <- condition$condition_id

  # Filter data for this condition
  cond_results <- results_df %>%
    filter(condition_id == cond_id) %>%
    mutate(
      fitted_model_code = factor(fitted_model_code,
        levels = c("NG", "NC", "SG", "SC")
      ),
      fit_type = ifelse(fitted_model_code == condition$correct_fit_model_code,
        "Correct", "Standard"
      )
    )

  cond_sampler <- sampler_info_df %>%
    filter(condition_id == cond_id) %>%
    mutate(
      fitted_model_code = factor(fitted_model_code,
        levels = c("NG", "NC", "SG", "SC")
      ),
      fit_type = ifelse(fitted_model_code == condition$correct_fit_model_code,
        "Correct", "Standard"
      )
    )

  # Skip if no data
  if (nrow(cond_results) == 0 && nrow(cond_sampler) == 0) {
    cat(sprintf("Skipping condition %03d: No data available.\n", cond_id))
    next
  }

  # Generate report
  generate_condition_report(
    condition_id = cond_id,
    condition_info = condition,
    condition_data = cond_results,
    sampler_data = cond_sampler,
    output_dir = PLOTS_DIR,
    fits_dir = FITS_DIR
  )
}

# Step 6: Create overall summary CSV
cat("\n--- Creating Overall Summary ---\n")

# Aggregate diagnostics across all conditions
overall_summary <- sampler_info_df %>%
  left_join(select(condition_info_df, condition_id, correct_fit_model_code),
    by = "condition_id"
  ) %>%
  mutate(fit_type = ifelse(fitted_model_code == correct_fit_model_code,
    "Correct", "Standard"
  )) %>%
  group_by(condition_id, fitted_model_code, fit_type) %>%
  summarise(
    n_reps = n(),
    pct_diverged = mean(divergences > 0, na.rm = TRUE) * 100,
    avg_divergences = mean(divergences, na.rm = TRUE),
    avg_efmi = mean(eFMI, na.rm = TRUE),
    pct_low_efmi = mean(eFMI < 0.3, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Save summary
write.csv(overall_summary,
  file.path(RESULTS_DIR, "convergence_summary_all_conditions.csv"),
  row.names = FALSE
)

cat("\n=== Convergence Analysis Complete ===\n")
cat(sprintf("Reports saved in: %s\n", PLOTS_DIR))
cat(sprintf("Overall summary saved in: %s\n", RESULTS_DIR))
