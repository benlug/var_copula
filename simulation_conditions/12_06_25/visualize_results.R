#!/usr/bin/env Rscript
###########################################################################
# visualize_results.R
#   • reads max_Rhat from CSV (fast)
#   • SD‑Bias correctly calculated per status group
###########################################################################

suppressPackageStartupMessages({
  library(dplyr);  library(tidyr);  library(readr)
  library(stringr); library(purrr); library(ggplot2)
  library(gridExtra); library(grid)
})

## ── paths ---------------------------------------------------------------
# Robust path setting
if (!requireNamespace("this.path", quietly = TRUE)) {
  message("Package 'this.path' not found. Assuming current working directory.")
  BASE <- getwd() 
} else {
  tryCatch({
    BASE <- this.path::this.dir()
  }, error = function(e) {
    message("this.path::this.dir() failed. Assuming current working directory.")
    BASE <<- getwd()
  })
}

DIR  <- list(
  data    = file.path(BASE, "data"),
  res     = file.path(BASE, "results"),
  plots   = file.path(BASE, "results", "plots"),
  plots_g = file.path(BASE, "results", "plots_global"))
dir.create(DIR$plots,   FALSE, TRUE)
dir.create(DIR$plots_g, FALSE, TRUE)

files <- list(
  cond   = file.path(DIR$res, "summary_conditions.csv"),
  rep    = file.path(DIR$res, "summary_replications.csv"),
  design = file.path(DIR$data, "sim_conditions_singlelevel.rds"))

# Check if input files exist
if (!all(file.exists(unlist(files)))) {
    message("Missing required input files (summaries or design RDS). Run analysis first.")
    # Exit if run as a script, otherwise return silently
    if (sys.nframe() == 0) q(status = 1) else return(invisible())
}


## ── data loading and preparation ----------------------------------------
design <- readRDS(files$design) |>
  select(condition_id, skew_level, direction, T, rho, VARset)

# Load condition-level summary
cond <- read_csv(files$cond, show_col_types = FALSE) |>
  left_join(design, by = "condition_id")

# Load replication-level summary
rep  <- read_csv(files$rep,  show_col_types = FALSE) |>
  # Filter out rows where fitting failed completely (no parameters)
  filter(!is.na(param)) |> 
  left_join(design, by = "condition_id")

if (nrow(rep) == 0) {
    message("Replication summary is empty or contains no valid parameters. Skipping visualization.")
    if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

## parameter order --------------------------------------------------------
# Define the order for visualization
param_levels <- c("omega[1]","omega[2]","alpha[1]","alpha[2]",
                  "sigma[1]","sigma[2]",
                  "mu[1]","mu[2]",
                  "phi11","phi12","phi21","phi22","rho")
# Apply factor levels for consistent plotting
rep$param  <- factor(rep$param,  levels = param_levels)
cond$param <- factor(cond$param, levels = param_levels)

## status (clean / problem) classification --------------------------------
# Classify runs based on MCMC diagnostics
RHAT_THRESHOLD <- 1.01
rep <- rep |>
  mutate(bad_rhat  = ifelse(is.na(max_rhat), FALSE, max_rhat > RHAT_THRESHOLD),
         divergent = n_div > 0,
         # "problem" if Rhat is high OR there are divergences
         status    = ifelse(bad_rhat | divergent, "problem", "clean"))

## plotting helpers -------------------------------------------------------
# Standard facet setup (Direction vs VAR set)
facet_basic  <- facet_grid(direction ~ VARset,
                           labeller = labeller(VARset = label_both))
# Facet setup including MCMC status
facet_status <- facet_grid(status + direction ~ VARset,
                           labeller = labeller(
                             VARset = label_both,
                             status = c(clean = "Clean Fits",
                                        problem = "Problematic Fits")))

# Helper to create table grobs for inclusion in PDFs
table_grob <- function(df) gridExtra::tableGrob(
  df, rows = NULL,
  theme = gridExtra::ttheme_minimal(base_size = 8,
                                    padding = unit(c(2,2), "mm")))

# Function to write a standard metric page (aggregated across all fits)
write_metric_page <- function(df, metric, ylab, pdf_file,
                              facet_fun, table = FALSE) {
  # Basic bar plot visualization
  g <- ggplot(df, aes(x = param, y = .data[[metric]], fill = model)) +
         geom_col(position = position_dodge(.6), width = .5) +
         coord_flip() + # Flip coordinates for better readability of params
         facet_fun + 
         theme_bw(8) +
         labs(x = "", y = ylab, title = paste("Global:", metric))
  
  # If the metric is coverage, add a reference line at 0.95
  if (metric == "coverage_95") {
      g <- g + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red")
  }
  
  pdf(pdf_file, 11, 7)
  # Optionally include a data table in the PDF
  if (table && metric %in% names(df)) {
    tbl_data <- select(df, condition_id, model, param, !!metric := .data[[metric]])
    gridExtra::grid.arrange(table_grob(tbl_data))
  }
  print(g)
  dev.off()
}

# Function to write metrics split by MCMC status (clean vs problem)
write_metric_split <- function(df, metric_col, ylab, pdf_file) {
  tmp <- df
  
  # Special handling for SD-Bias: needs recalculation within status groups
  if (metric_col == "sd_bias") {
    # Calculate Empirical SD and SD-Bias *within* each status group
    tmp <- tmp |>
      group_by(condition_id, model, param, status) |>
      mutate(
          N_group = n(),
          # Ensure sufficient samples for SD calculation (sd needs >1 sample)
          # If N=1, empirical SD is treated as 0 for bias calculation.
          emp_sd = ifelse(N_group > 1, sd(post_mean), 0), 
          .val   = post_sd - emp_sd
      ) |>
      ungroup()
  } else {
    # For other metrics, just use the column value directly
    if (metric_col %in% names(tmp)) {
        tmp$.val <- tmp[[metric_col]]
    } else {
        return(invisible()) # Skip if metric column doesn't exist
    }
  }

  # Count number of reps in each status group (for context)
  cnt <- tmp |> count(status, direction, VARset, name = "n_reps")
  
  # Aggregate the metric value (mean) within each group
  agg <- tmp |>
    group_by(condition_id, model, param, status, direction, VARset) |>
    summarise(mean_val = mean(.data$.val, na.rm = TRUE), .groups = "drop")

  # Visualization
  g <- ggplot(agg, aes(x = param, y = mean_val, fill = model)) +
      geom_col(position = "dodge", width = .5) +
      coord_flip() + 
      facet_status + 
      theme_bw(8) +
      labs(y = ylab, x = "", title = paste("Split by Status:", metric_col))
      
  # Add reference line for coverage
  if (metric_col == "cover95") {
      g <- g + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red")
  }

  pdf(pdf_file, 11, 8.5) # slightly taller for split plots
  # Include the counts table
  gridExtra::grid.arrange(table_grob(cnt))
  print(g)
  dev.off()
}

# Function to summarize divergences
write_div_table <- function(df, pdf_file) {
  # Summarize total divergences and count reps with divergences
  t <- df |> group_by(condition_id, model, direction, VARset, T, rho) |>
           summarise(Total_Divergences = sum(n_div), 
                     N_Reps = n(), 
                     Reps_with_Div = sum(n_div > 0),
                     .groups = "drop")
  pdf(pdf_file, 11, 7)
  gridExtra::grid.arrange(table_grob(t))
  dev.off()
}

# Function for global overview plots (across T, faceted by other factors)
global_plots <- function(df, metric, ylab, sk_tag) {
  if (!metric %in% names(df)) return(invisible())

  # Iterate over rho values
  for (r in sort(unique(df$rho))) {
    d <- filter(df, rho == r)
    tag <- paste0(sk_tag, "_rho", r, "_", metric)
    
    # Line plot showing trends over T
    g1 <- ggplot(d, aes(x = factor(T), y = .data[[metric]], colour = model,
                        group  = model)) +
            geom_line() + geom_point() +
            # Facet by parameter, direction, and VARset; allow free Y scales
            facet_grid(param ~ direction + VARset, scales = "free_y",
                       labeller = labeller(VARset = label_both)) +
            theme_bw(7) +
            labs(x = "Time Points (T)", y = ylab, 
                 title = paste0(metric, " Trends | ρ=", r))
    
    if (metric == "coverage_95") {
        g1 <- g1 + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red")
    }

    pdf(file.path(DIR$plots_g, paste0(tag, ".pdf")), 12, 9)
    print(g1)
    dev.off()
  }
}

## glossary ---------------------------------------------------------------
# Write definitions to a text file for reference
writeLines(c(
  "Metric definitions:",
  "  relative bias  = (posterior mean − truth) / |truth| (if truth!=0, else bias)",
  "  coverage_95    = proportion of reps where 95% CI covers truth",
  "  posterior SD   = mean posterior SD across reps",
  "  SD‑Bias        = Mean(Posterior SD) − Empirical SD(Posterior Means)",
  "  n_div          = divergent transitions (post‑warm‑up)",
  paste("  clean          = no divergence & Rhat ≤", RHAT_THRESHOLD)),
  file.path(DIR$res, "metric_glossary.txt"))

## metric map -------------------------------------------------------------
# Mapping between condition-level (aggregated) names and replication-level names
metrics <- list(
  # Aggregated Name    # Replication Name   # Y-axis Label
  mean_rel_bias = c("rel_bias",  "Relative Bias"),
  coverage_95   = c("cover95",   "95% Coverage Proportion"),
  mean_post_sd  = c("post_sd",   "Mean Posterior SD"),
  sd_bias       = c("sd_bias",   "SD-Bias"),
  mean_n_div    = c("n_div",     "Mean # Divergences"))

## main loop --------------------------------------------------------------
message("Starting visualization...")
# Iterate through each skew level
for (sk in unique(design$skew_level)) {
  message("Processing skew level: ", sk)
  sk_dir <- file.path(DIR$plots, sk); dir.create(sk_dir, FALSE, TRUE)
  
  # Filter data for the current skew level
  cond_sk <- filter(cond, skew_level == sk)
  rep_sk  <- filter(rep,  skew_level == sk)

  # Iterate through unique combinations of T and rho
  combos <- unique(cond_sk[c("T","rho")])
  for (i in seq_len(nrow(combos))) {
    Tval <- combos$T[i]; rho_val <- combos$rho[i]
    tag  <- paste0("T", Tval, "_rho", rho_val)
    outd <- file.path(sk_dir, tag); dir.create(outd, FALSE, TRUE)

    # Further filter data for the specific T/rho combination
    sel_c <- filter(cond_sk, T == Tval, rho == rho_val)
    sel_r <- filter(rep_sk,  T == Tval, rho == rho_val)

    # Generate plots for each metric
    walk(names(metrics), function(mc_agg) {
      mc_rep <- metrics[[mc_agg]][1]
      ylab <- metrics[[mc_agg]][2]
      
      # 1. Global plots (all fits aggregated)
      write_metric_page(
        sel_c, mc_agg, ylab,
        file.path(outd, paste0(mc_agg, ".pdf")),
        facet_basic, table = TRUE)
        
      # 2. Split plots (by clean/problem status)
      write_metric_split(
        sel_r, mc_rep, ylab,
        file.path(outd, paste0(mc_agg, "_split.pdf")))
    })

    # 3. Divergence summary table
    write_div_table(sel_r, file.path(outd, "divergence_table.pdf"))
  }

  # Generate global overview plots for this skew level
  walk(names(metrics), function(mc) {
    global_plots(cond_sk, mc, metrics[[mc]][2], sk)
  })
}

message("✓ visualisation complete – see ", DIR$plots, " & ", DIR$plots_g)