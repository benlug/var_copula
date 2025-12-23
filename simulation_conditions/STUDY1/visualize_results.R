#!/usr/bin/env Rscript
###########################################################################
# visualize_results.R
# Improved visualization with cleaner aesthetics and organization
# FIXED: n_reps now counts distinct replications, not parameter rows
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

## ── Paths ----------------------------------------------------------------
if (!requireNamespace("this.path", quietly = TRUE)) {
  BASE <- getwd()
} else {
  tryCatch({
    BASE <- this.path::this.dir()
  }, error = function(e) {
    BASE <<- getwd()
  })
}

DIR <- list(
  data   = file.path(BASE, "data"),
  res    = file.path(BASE, "results"),
  plots  = file.path(BASE, "results", "plots")
)
dir.create(DIR$plots, FALSE, TRUE)

files <- list(
  cond   = file.path(DIR$res, "summary_conditions.csv"),
  rep    = file.path(DIR$res, "summary_replications.csv"),
  design = file.path(DIR$data, "sim_conditions_singlelevel.rds")
)

if (!all(file.exists(unlist(files)))) {
  message("Missing required input files. Run analysis first.")
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

## ── Data Loading ---------------------------------------------------------
design <- readRDS(files$design) |>
  select(condition_id, skew_level, direction, T, rho, VARset)

cond <- read_csv(files$cond, show_col_types = FALSE) |>
  left_join(design, by = "condition_id")

rep <- read_csv(files$rep, show_col_types = FALSE) |>
  filter(!is.na(param)) |>
  left_join(design, by = "condition_id")

if (nrow(rep) == 0) {
  message("No valid data. Skipping visualization.")
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

## ── Factor Setup ---------------------------------------------------------
param_levels <- c("mu[1]", "mu[2]",
                  "phi11", "phi12", "phi21", "phi22",
                  "rho",
                  "sigma[1]", "sigma[2]",
                  "omega[1]", "omega[2]", "alpha[1]", "alpha[2]")

rep$param  <- factor(rep$param,  levels = param_levels)
cond$param <- factor(cond$param, levels = param_levels)

rep$model  <- factor(rep$model,  levels = c("NG", "SG"))
cond$model <- factor(cond$model, levels = c("NG", "SG"))

rep$skew_level  <- factor(rep$skew_level,  levels = c("moderateSN", "strongSN", "extremeCHI"))
cond$skew_level <- factor(cond$skew_level, levels = c("moderateSN", "strongSN", "extremeCHI"))

## ── MCMC Status Classification -------------------------------------------
RHAT_THRESHOLD <- 1.01
rep <- rep |>
  mutate(
    mcmc_status = case_when(
      is.na(max_rhat) | status != "ok" ~ "Failed",
      max_rhat > RHAT_THRESHOLD | n_div > 0 ~ "Problematic",
      TRUE ~ "Clean"
    ),
    mcmc_status = factor(mcmc_status, levels = c("Clean", "Problematic", "Failed"))
  )

## ── Theme Setup ----------------------------------------------------------
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", size = rel(0.85)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(color = "grey40", size = rel(0.9)),
      axis.title = element_text(size = rel(0.95))
    )
}

# Color palette
model_colors <- c("NG" = "#E69F00", "SG" = "#56B4E9")
status_colors <- c("Clean" = "#4DAF4A", "Problematic" = "#FF7F00", "Failed" = "#E41A1C")

## ── Core Parameters ------------------------------------------------------
core_params <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")

## ── Helper Functions -----------------------------------------------------

# Calculate RMSE
cond <- cond |>
  mutate(RMSE = sqrt(mean_bias^2 + coalesce(emp_sd^2, 0)))

# FIXED: Count distinct replications, not parameter rows
count_replications <- function(df) {
  df |>
    distinct(condition_id, rep_id, model, mcmc_status, T, skew_level, direction, VARset, rho) |>
    group_by(model, T, skew_level, direction, VARset, rho, mcmc_status) |>
    summarise(n_reps = n(), .groups = "drop")
}

# Aggregate by MCMC status for split analyses
aggregate_by_status <- function(df) {
  df |>
    filter(mcmc_status != "Failed") |>
    group_by(condition_id, model, param, mcmc_status, T, skew_level, direction, VARset, rho) |>
    summarise(
      N_valid = n(),
      mean_rel_bias = mean(rel_bias, na.rm = TRUE),
      coverage_95 = mean(cover95, na.rm = TRUE),
      mean_post_sd = mean(post_sd, na.rm = TRUE),
      emp_sd = if(n() > 1) sd(post_mean, na.rm = TRUE) else 0,
      mean_bias = mean(bias, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      sd_bias = mean_post_sd - emp_sd,
      RMSE = sqrt(mean_bias^2 + emp_sd^2)
    )
}

## ── Plot Functions -------------------------------------------------------

# 1. Bias plot by skewness condition
plot_bias_by_condition <- function(data, skew_lvl) {
  d <- data |>
    filter(skew_level == skew_lvl, param %in% core_params)
  
  ggplot(d, aes(x = factor(T), y = mean_rel_bias, color = model, group = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_grid(param ~ direction + VARset, 
               scales = "free_y",
               labeller = labeller(VARset = label_both)) +
    scale_color_manual(values = model_colors) +
    theme_publication() +
    labs(
      title = sprintf("Relative Bias: %s", skew_lvl),
      subtitle = "Dashed line = zero bias. Note: mu panels show absolute bias (truth=0)",
      x = "Time Points (T)",
      y = "Mean Relative Bias",
      color = "Model"
    )
}

# 2. Coverage plot by skewness condition
plot_coverage_by_condition <- function(data, skew_lvl) {
  d <- data |>
    filter(skew_level == skew_lvl, param %in% core_params)
  
  ylims <- if (skew_lvl == "extremeCHI") c(0.5, 1.0) else c(0.8, 1.0)
  
  ggplot(d, aes(x = factor(T), y = coverage_95, color = model, group = model)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_grid(param ~ direction + VARset,
               labeller = labeller(VARset = label_both)) +
    scale_color_manual(values = model_colors) +
    coord_cartesian(ylim = ylims) +
    theme_publication() +
    labs(
      title = sprintf("95%% Coverage: %s", skew_lvl),
      subtitle = "Dashed line = nominal 0.95 coverage",
      x = "Time Points (T)",
      y = "Empirical Coverage",
      color = "Model"
    )
}

# 3. SD-Bias plot
plot_sdbias_by_condition <- function(data, skew_lvl) {
  d <- data |>
    filter(skew_level == skew_lvl, param %in% core_params)
  
  ggplot(d, aes(x = factor(T), y = sd_bias, color = model, group = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_grid(param ~ direction + VARset,
               scales = "free_y",
               labeller = labeller(VARset = label_both)) +
    scale_color_manual(values = model_colors) +
    theme_publication() +
    labs(
      title = sprintf("SD-Bias: %s", skew_lvl),
      subtitle = "Mean(Posterior SD) - Empirical SD | Dashed line = 0",
      x = "Time Points (T)",
      y = "SD-Bias",
      color = "Model"
    )
}

# 4. MCMC status overview - FIXED to count replications correctly
plot_mcmc_status <- function(rep_data) {
  # Count distinct replications, not parameter rows
  summary_data <- rep_data |>
    distinct(condition_id, rep_id, model, mcmc_status, T, skew_level) |>
    group_by(model, T, skew_level, mcmc_status) |>
    summarise(Count = n(), .groups = "drop")
  
  ggplot(summary_data, aes(x = factor(T), y = Count, fill = mcmc_status)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(model ~ skew_level) +
    scale_fill_manual(values = status_colors) +
    theme_publication() +
    labs(
      title = "MCMC Convergence Status",
      subtitle = sprintf("Classification: Clean (R-hat ≤ %.2f, no divergences) vs Problematic. Counts are distinct replications.", RHAT_THRESHOLD),
      x = "Time Points (T)",
      y = "Number of Replications",
      fill = "Status"
    )
}

# 5. Divergence distribution
plot_divergence_dist <- function(rep_data) {
  # Get one row per replication (use first parameter, e.g., rho)
  div_data <- rep_data |>
    filter(param == "rho", mcmc_status != "Failed") |>
    distinct(condition_id, rep_id, model, T, skew_level, n_div)
  
  ggplot(div_data, aes(x = factor(T), y = n_div, fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
    facet_grid(model ~ skew_level) +
    scale_fill_manual(values = model_colors) +
    theme_publication() +
    labs(
      title = "Distribution of Divergent Transitions",
      subtitle = "Post-warmup divergences per replication",
      x = "Time Points (T)",
      y = "Divergence Count",
      fill = "Model"
    )
}

# 6. Coverage by MCMC status - FIXED counting
plot_coverage_by_status <- function(rep_data, model_name) {
  status_agg <- aggregate_by_status(rep_data)
  
  d <- status_agg |>
    filter(model == model_name, param %in% core_params)
  
  overview <- d |>
    group_by(T, param, skew_level, mcmc_status) |>
    summarise(mean_coverage = mean(coverage_95, na.rm = TRUE), .groups = "drop")
  
  # Add count annotation
  counts <- count_replications(rep_data) |>
    filter(model == model_name) |>
    group_by(T, skew_level, mcmc_status) |>
    summarise(total_reps = sum(n_reps), .groups = "drop")
  
  ggplot(overview, aes(x = factor(T), y = mean_coverage, 
                       color = mcmc_status, group = mcmc_status)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_grid(param ~ skew_level) +
    scale_color_manual(values = status_colors[1:2]) +
    coord_cartesian(ylim = c(0.75, 1.0)) +
    theme_publication() +
    labs(
      title = sprintf("%s Model: Coverage by MCMC Status", model_name),
      subtitle = "Averaged over VAR set, rho, and direction",
      x = "Time Points (T)",
      y = "Mean Coverage",
      color = "MCMC Status"
    )
}

# 7. Marginal parameter analysis
plot_marginal_params <- function(cond_data, model_name) {
  if (model_name == "NG") {
    params <- c("sigma[1]", "sigma[2]")
    title <- "NG Model: Scale Parameter (σ) Bias"
    subtitle <- "True value = 1 (standardized innovations)"
  } else {
    params <- c("alpha[1]", "alpha[2]")
    title <- "SG Model: Shape Parameter (α) Relative Bias"
    subtitle <- "Only defined for skew-normal DGP conditions"
  }
  
  d <- cond_data |>
    filter(model == model_name, param %in% params)
  
  if (model_name == "SG") {
    d <- d |> filter(skew_level %in% c("moderateSN", "strongSN"))
  }
  
  if (nrow(d) == 0) return(NULL)
  
  ggplot(d, aes(x = factor(T), y = mean_rel_bias, color = skew_level, group = skew_level)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_grid(param ~ direction + VARset, 
               scales = "free_y",
               labeller = labeller(VARset = label_both)) +
    theme_publication() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time Points (T)",
      y = "Relative Bias",
      color = "DGP"
    )
}

# 8. Replication counts summary table - NEW
create_rep_count_summary <- function(rep_data) {
  rep_data |>
    distinct(condition_id, rep_id, model, mcmc_status, T, skew_level, direction, VARset, rho) |>
    group_by(model, skew_level, T, mcmc_status) |>
    summarise(n_replications = n(), .groups = "drop") |>
    pivot_wider(names_from = mcmc_status, values_from = n_replications, values_fill = 0) |>
    mutate(Total = Clean + Problematic + Failed)
}

## ── Generate All Plots ---------------------------------------------------
message("Generating plots...")

# Create subdirectories
for (sk in c("moderateSN", "strongSN", "extremeCHI")) {
  dir.create(file.path(DIR$plots, sk), FALSE, TRUE)
}
dir.create(file.path(DIR$plots, "overview"), FALSE, TRUE)

# Save replication count summary
rep_counts <- create_rep_count_summary(rep)
write_csv(rep_counts, file.path(DIR$res, "replication_counts_by_status.csv"))
message("  - Replication counts saved to replication_counts_by_status.csv")

# Overview plots
ggsave(file.path(DIR$plots, "overview", "mcmc_status.pdf"),
       plot_mcmc_status(rep), width = 12, height = 8)
message("  - MCMC status overview")

ggsave(file.path(DIR$plots, "overview", "divergence_distribution.pdf"),
       plot_divergence_dist(rep), width = 12, height = 8)
message("  - Divergence distribution")

ggsave(file.path(DIR$plots, "overview", "coverage_by_status_SG.pdf"),
       plot_coverage_by_status(rep, "SG"), width = 12, height = 10)
message("  - SG coverage by status")

ggsave(file.path(DIR$plots, "overview", "coverage_by_status_NG.pdf"),
       plot_coverage_by_status(rep, "NG"), width = 12, height = 10)
message("  - NG coverage by status")

ggsave(file.path(DIR$plots, "overview", "marginal_sigma_NG.pdf"),
       plot_marginal_params(cond, "NG"), width = 14, height = 6)
message("  - NG sigma parameters")

p <- plot_marginal_params(cond, "SG")
if (!is.null(p)) {
  ggsave(file.path(DIR$plots, "overview", "marginal_alpha_SG.pdf"), p, width = 14, height = 6)
  message("  - SG alpha parameters")
}

# Per-condition plots
for (sk in c("moderateSN", "strongSN", "extremeCHI")) {
  message(sprintf("  - %s plots...", sk))
  
  ggsave(file.path(DIR$plots, sk, "bias.pdf"),
         plot_bias_by_condition(cond, sk), width = 16, height = 12)
  
  ggsave(file.path(DIR$plots, sk, "coverage.pdf"),
         plot_coverage_by_condition(cond, sk), width = 16, height = 12)
  
  ggsave(file.path(DIR$plots, sk, "sd_bias.pdf"),
         plot_sdbias_by_condition(cond, sk), width = 16, height = 12)
}

message("\n=== Visualization Summary ===")
message("Total distinct replications analyzed:")
print(rep_counts |> group_by(model, skew_level) |> summarise(Total = sum(Total), .groups = "drop"))

message("\nVisualization complete. Plots saved to: ", DIR$plots)
