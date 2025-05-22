#!/usr/bin/env Rscript
# Identify which simulation replications failed (missing data or fit files)

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
# Use optional base directory argument; default current working directory
base_dir <- if (length(args) >= 1) args[1] else getwd()

# Paths follow the directory layout from run_pipeline.R
DATA_DIR <- file.path(base_dir, "../data")
FITS_DIR <- file.path(base_dir, "../fits")
conditions_file <- file.path(DATA_DIR, "sim_conditions.rds")

if (!file.exists(conditions_file)) {
  stop("Simulation conditions file not found: ", conditions_file)
}

sim_conditions <- readRDS(conditions_file)

missing_records <- list()

for (i in seq_len(nrow(sim_conditions))) {
  cond <- sim_conditions[i, ]
  cond_id <- as.integer(cond$condition_id)
  n_reps <- as.integer(cond$n_reps)
  missing_reps <- integer(0)
  for (r in seq_len(n_reps)) {
    data_file <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cond_id, r))
    if (!file.exists(data_file)) {
      missing_reps <- c(missing_reps, r)
    }
  }
  if (length(missing_reps) > 0) {
    missing_records[[as.character(cond_id)]] <- missing_reps
  }
}

if (length(missing_records) == 0) {
  cat("All expected replication files are present.\n")
} else {
  cat("Missing replication files by condition:\n")
  for (cid in names(missing_records)) {
    reps <- paste(missing_records[[cid]], collapse = ", ")
    cat(sprintf("  Condition %s -> Replications: %s\n", cid, reps))
  }
}
