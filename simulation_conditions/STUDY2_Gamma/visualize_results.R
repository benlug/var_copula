#!/usr/bin/env Rscript
###########################################################################
# visualize_results.R (Updated for Study 1 and Study 2 compatibility)
#   • Dynamic directory detection (S1 vs S2)
#   • Handles EG parameters
#   • Robust to different design grid column names
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

## ── paths ---------------------------------------------------------------
# Set BASE directory robustly
if (!requireNamespace("this.path", quietly = TRUE)) {
  BASE <- getwd()
} else {
  tryCatch(
    {
      # Try to find the directory of the script if sourced/run via Rscript
      BASE <- this.path::this.dir()
    },
    error = function(e) {
      # Fallback to current working directory if this.path fails
      BASE <<- getwd()
    }
  )
}

# Dynamic Directory Detection (Study 1 vs Study 2)
# Checks for the existence of Study 2 results directory
if (dir.exists(file.path(BASE, "results"))) {
  STUDY_TAG <- "S2"
  message("Visualizing Study 2 results.")
  DIR <- list(
    data    = file.path(BASE, "data"),
    res     = file.path(BASE, "results"),
    plots   = file.path(BASE, "results", "plots"),
    plots_g = file.path(BASE, "results", "plots_global")
  )
  DESIGN_FILE <- "sim_conditions.rds"
} else {
  STUDY_TAG <- "S1"
  message("Visualizing Study 1 results.")
  # Fallback to Study 1 directory structure
  DIR <- list(
    data    = file.path(BASE, "data"),
    res     = file.path(BASE, "results"),
    plots   = file.path(BASE, "results", "plots"),
    plots_g = file.path(BASE, "results", "plots_global")
  )
  DESIGN_FILE <- "sim_conditions_singlelevel.rds"
}

dir.create(DIR$plots, FALSE, TRUE)
dir.create(DIR$plots_g, FALSE, TRUE)

files <- list(
  cond   = file.path(DIR$res, "summary_conditions.csv"),
  rep    = file.path(DIR$res, "summary_replications.csv"),
  design = file.path(DIR$data, DESIGN_FILE)
)

# Check if input files exist
if (!all(file.exists(unlist(files)))) {
  message("Missing required input files (summaries or design RDS). Run analysis first.")
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}


## ── data loading and preparation ----------------------------------------
design <- readRDS(files$design)

# Harmonize column names between S1 (skew_level) and S2 (dgp_level)
# We standardize on 'skew_level' internally for consistency.
if ("dgp_level" %in% names(design) && !"skew_level" %in% names(design)) {
  names(design)[names(design) == "dgp_level"] <- "skew_level"
}

# Select relevant columns, ensuring they exist
required_cols <- c("condition_id", "skew_level", "direction", "T", "rho", "VARset")
if (!all(required_cols %in% names(design))) {
  stop("Design file is missing required columns.")
}
design <- select(design, all_of(required_cols))


# Load condition-level summary
cond <- read_csv(files$cond, show_col_types = FALSE) |>
  left_join(design, by = "condition_id")

# Load replication-level summary
rep <- read_csv(files$rep, show_col_types = FALSE) |>
  # Filter out rows where fitting failed completely (no parameters)
  filter(!is.na(param)) |>
  left_join(design, by = "condition_id")

if (nrow(rep) == 0) {
  message("Replication summary is empty or contains no valid parameters. Skipping visualization.")
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

## parameter order --------------------------------------------------------
# Define the order for visualization (Updated for Study 2)
param_levels <- c(
  # Study 2 (EG)
  "sigma_exp[1]", "sigma_exp[2]",
  # Study 3 (GG)
  "sigma_gam[1]", "sigma_gam[2]",
  # Study 1 (SG)
  "omega[1]", "omega[2]", "alpha[1]", "alpha[2]",
  # Study 1 (NG)
  "sigma[1]", "sigma[2]",
  # Core parameters
  "mu[1]", "mu[2]",
  "phi11", "phi12", "phi21", "phi22", "rho"
)

# Apply factor levels for consistent plotting.
# Filter parameters that actually exist in the data to avoid empty levels in plots.
# This handles cases where Study 1 runs won't have sigma_exp, etc.
existing_params <- intersect(param_levels, unique(c(levels(rep$param), unique(rep$param))))

rep$param <- factor(rep$param, levels = existing_params)
cond$param <- factor(cond$param, levels = existing_params)

## status (clean / problem) classification --------------------------------
# Classify runs based on MCMC diagnostics
RHAT_THRESHOLD <- 1.01
rep <- rep |>
  mutate(
    bad_rhat = ifelse(is.na(max_rhat), FALSE, max_rhat > RHAT_THRESHOLD),
    # Ensure n_div is treated as 0 if NA
    n_div_clean = ifelse(is.na(n_div), 0, n_div),
    divergent = n_div_clean > 0,
    # "problem" if Rhat is high OR there are divergences
    status = ifelse(bad_rhat | divergent, "problem", "clean")
  )

## plotting helpers -------------------------------------------------------

# Standard facet setup (Direction vs VAR set)
facet_basic <- facet_grid(direction ~ VARset,
  labeller = labeller(VARset = label_both)
)
# Facet setup including MCMC status
facet_status <- facet_grid(status + direction ~ VARset,
  labeller = labeller(
    VARset = label_both,
    status = c(
      clean = "Clean Fits",
      problem = "Problematic Fits"
    )
  )
)

# Helper to create table grobs for inclusion in PDFs
table_grob <- function(df) {
  gridExtra::tableGrob(
    df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = unit(c(2, 2), "mm")
    )
  )
}

# Function to write a standard metric page (aggregated across all fits)
write_metric_page <- function(df, metric, ylab, pdf_file,
                              facet_fun, table = FALSE) {
  # Filter out NAs for the metric (e.g., bias when truth is NA due to misspecification)
  df_plot <- filter(df, !is.na(.data[[metric]]))
  if (nrow(df_plot) == 0) {
    return(invisible())
  }

  # Basic bar plot visualization
  g <- ggplot(df_plot, aes(x = param, y = .data[[metric]], fill = model)) +
    geom_col(position = position_dodge(.6), width = .5) +
    coord_flip() +
    facet_fun +
    theme_bw(8) +
    labs(x = "", y = ylab, title = paste("Global:", metric))

  # If the metric is coverage, add a reference line at 0.95
  if (metric == "coverage_95") {
    g <- g + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red")
  }

  # Use tryCatch for robust PDF generation
  tryCatch({
    pdf(pdf_file, 11, 7)
    # Optionally include a data table in the PDF
    if (table && metric %in% names(df) && nrow(df_plot) > 0) {
      # Select relevant columns and pivot wider for easier reading in the table
      tbl_data <- df_plot |>
        select(condition_id, model, param, value = .data[[metric]]) |>
        mutate(value = round(value, 3)) |>
        # Requires tidyr >= 1.0.0
        pivot_wider(names_from = model, values_from = value) |>
        arrange(condition_id, param)

      # Split table if it's too long (PDF Pagination)
      n_rows <- nrow(tbl_data)
      rows_per_page <- 40
      n_pages <- ceiling(n_rows / rows_per_page)

      for (p in 1:n_pages) {
        start_row <- (p - 1) * rows_per_page + 1
        end_row <- min(p * rows_per_page, n_rows)
        gridExtra::grid.arrange(table_grob(tbl_data[start_row:end_row, ]))
      }
    }
    print(g)
  }, error = function(e) {
    message("Error writing PDF ", pdf_file, ": ", e$message)
  }, finally = {
    if (names(dev.cur()) != "null device") dev.off()
  })
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
        # If N=1, empirical SD is treated as 0 for bias calculation.
        emp_sd = ifelse(N_group > 1, sd(post_mean, na.rm = TRUE), 0),
        .val = post_sd - emp_sd
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

  # Filter out NAs (e.g., bias/coverage when truth is NA)
  tmp <- filter(tmp, !is.na(.val))
  if (nrow(tmp) == 0) {
    return(invisible())
  }

  # Count number of reps in each status group (for context)
  cnt <- tmp |>
    count(status, direction, VARset, name = "n_reps") |>
    pivot_wider(names_from = status, values_from = n_reps, values_fill = list(n_reps = 0))

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

  tryCatch({
    pdf(pdf_file, 11, 8.5) # slightly taller for split plots
    # Include the counts table
    gridExtra::grid.arrange(table_grob(cnt))
    print(g)
  }, error = function(e) {
    message("Error writing PDF ", pdf_file, ": ", e$message)
  }, finally = {
    if (names(dev.cur()) != "null device") dev.off()
  })
}

# Function to summarize divergences
write_div_table <- function(df, pdf_file) {
  # Summarize total divergences and count reps with divergences
  t <- df |>
    group_by(condition_id, model, direction, VARset, T, rho) |>
    # Use the cleaned n_div count
    summarise(
      Total_Divergences = sum(n_div_clean, na.rm = TRUE),
      N_Reps = n(),
      Reps_with_Div = sum(divergent, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(Total_Divergences))

  tryCatch({
    pdf(pdf_file, 11, 7)
    # Split table if too long (PDF Pagination)
    n_rows <- nrow(t)
    rows_per_page <- 50
    n_pages <- ceiling(n_rows / rows_per_page)

    if (n_pages > 0) {
      for (p in 1:n_pages) {
        start_row <- (p - 1) * rows_per_page + 1
        end_row <- min(p * rows_per_page, n_rows)
        gridExtra::grid.arrange(table_grob(t[start_row:end_row, ]))
      }
    }
  }, error = function(e) {
    message("Error writing PDF ", pdf_file, ": ", e$message)
  }, finally = {
    if (names(dev.cur()) != "null device") dev.off()
  })
}

# Function for global overview plots (across T, faceted by other factors)
global_plots <- function(df, metric, ylab, sk_tag) {
  if (!metric %in% names(df)) {
    return(invisible())
  }

  # Filter out NAs
  df_plot <- filter(df, !is.na(.data[[metric]]))
  if (nrow(df_plot) == 0) {
    return(invisible())
  }

  # Iterate over rho values
  for (r in sort(unique(df_plot$rho))) {
    d <- filter(df_plot, rho == r)
    tag <- paste0(sk_tag, "_rho", r, "_", metric)

    # Line plot showing trends over T
    g1 <- ggplot(d, aes(
      x = factor(T), y = .data[[metric]], colour = model,
      group = model
    )) +
      geom_line() +
      geom_point() +
      # Facet by parameter, direction, and VARset; allow free Y scales
      # Use drop=TRUE to remove empty facets if parameters don't exist for all models
      facet_grid(param ~ direction + VARset,
        scales = "free_y", drop = TRUE,
        labeller = labeller(VARset = label_both)
      ) +
      theme_bw(7) +
      labs(
        x = "Time Points (T)", y = ylab,
        title = paste0(metric, " Trends | ρ=", r)
      )

    if (metric == "coverage_95") {
      g1 <- g1 + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red")
    }

    pdf_file <- file.path(DIR$plots_g, paste0(tag, ".pdf"))
    tryCatch({
      # Adjust height based on the number of parameters shown
      n_params <- length(unique(d$param))
      # Dynamic height calculation
      pdf_height <- max(5, min(15, n_params * 0.8))

      pdf(pdf_file, 12, pdf_height)
      print(g1)
    }, error = function(e) {
      message("Error writing PDF ", pdf_file, ": ", e$message)
    }, finally = {
      if (names(dev.cur()) != "null device") dev.off()
    })
  }
}

## glossary ---------------------------------------------------------------
# Write definitions to a text file for reference
writeLines(
  c(
    paste("Study:", STUDY_TAG),
    "Metric definitions:",
    "  relative bias  = (posterior mean − truth) / |truth| (if truth!=0, else bias). NA if truth unavailable.",
    "  coverage_95    = proportion of reps where 95% CI covers truth. NA if truth unavailable.",
    "  posterior SD   = mean posterior SD across reps",
    "  SD‑Bias        = Mean(Posterior SD) − Empirical SD(Posterior Means)",
    "  n_div          = divergent transitions (post‑warm‑up)",
    paste("  clean          = no divergence & Rhat ≤", RHAT_THRESHOLD)
  ),
  file.path(DIR$res, "metric_glossary.txt")
)

## metric map -------------------------------------------------------------
# Mapping between condition-level (aggregated) names and replication-level names
metrics <- list(
  # Aggregated Name    # Replication Name   # Y-axis Label
  mean_rel_bias = c("rel_bias", "Relative Bias"),
  coverage_95   = c("cover95", "95% Coverage Proportion"),
  mean_post_sd  = c("post_sd", "Mean Posterior SD"),
  sd_bias       = c("sd_bias", "SD-Bias"),
  mean_n_div    = c("n_div", "Mean # Divergences")
)

## main loop --------------------------------------------------------------
message("Starting visualization...")
# Iterate through each skew level (or DGP level in S2)
# Use make.names to ensure the level is safe for directory names
design$skew_level_safe <- make.names(design$skew_level)

for (sk_safe in unique(design$skew_level_safe)) {
  # Get the original name for titles/labels
  sk_original <- unique(filter(design, skew_level_safe == sk_safe)$skew_level)

  message("Processing level: ", sk_original)
  sk_dir <- file.path(DIR$plots, sk_safe)
  dir.create(sk_dir, FALSE, TRUE)

  # Filter data for the current skew level
  cond_sk <- filter(cond, skew_level == sk_original)
  rep_sk <- filter(rep, skew_level == sk_original)

  if (nrow(cond_sk) == 0 || nrow(rep_sk) == 0) {
    message("  Skipping level ", sk_original, " (no data).")
    next
  }

  # Iterate through unique combinations of T and rho
  combos <- unique(cond_sk[c("T", "rho")])
  for (i in seq_len(nrow(combos))) {
    Tval <- combos$T[i]
    rho_val <- combos$rho[i]
    tag <- paste0("T", Tval, "_rho", rho_val)
    outd <- file.path(sk_dir, tag)
    dir.create(outd, FALSE, TRUE)

    # Further filter data for the specific T/rho combination
    sel_c <- filter(cond_sk, T == Tval, rho == rho_val)
    sel_r <- filter(rep_sk, T == Tval, rho == rho_val)

    if (nrow(sel_c) == 0 || nrow(sel_r) == 0) next

    # Generate plots for each metric
    walk(names(metrics), function(mc_agg) {
      mc_rep <- metrics[[mc_agg]][1]
      ylab <- metrics[[mc_agg]][2]

      # 1. Global plots (all fits aggregated)
      write_metric_page(
        sel_c, mc_agg, ylab,
        file.path(outd, paste0(mc_agg, ".pdf")),
        facet_basic,
        table = TRUE
      )

      # 2. Split plots (by clean/problem status)
      write_metric_split(
        sel_r, mc_rep, ylab,
        file.path(outd, paste0(mc_agg, "_split.pdf"))
      )
    })

    # 3. Divergence summary table
    write_div_table(sel_r, file.path(outd, "divergence_table.pdf"))
  }

  # Generate global overview plots for this skew level
  walk(names(metrics), function(mc) {
    # Pass the safe name for the filename
    global_plots(cond_sk, mc, metrics[[mc]][2], sk_safe)
  })
}

message("✓ visualisation complete – see ", DIR$plots, " & ", DIR$plots_g)
