#!/usr/bin/env Rscript
###########################################################################
# run_pipeline_gamma.R  – Study 2 (Gamma extension)
#
# Goal
#   Replicate Study 2's factorial design, replacing Exponential innovations
#   with standardized Gamma innovations.
#
# Update (Skewness levels)
#   This pipeline supports THREE Gamma-shape levels to vary skewness magnitude.
#   Skewness for Gamma(k, ·) is 2/sqrt(k); larger k = less skew.
#   The default set keeps the original k=2 as the LEAST-skewed level.
#
# Standardization
#   X ~ Gamma(shape=k, rate=beta)
#   Z = (X - k/beta) / (sqrt(k)/beta)
#
# Implementation note
#   We use beta = sqrt(k) by default for each k (so Var(X)=1), and then
#   Z = X - sqrt(k). Mirroring ("-") multiplies Z by -1 to induce left skew.
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  if (!requireNamespace("this.path", quietly = TRUE)) {
    if (interactive()) {
      install.packages("this.path")
    }
    library(this.path)
  }
})

## -- folders -------------------------------------------------------------
tryCatch(
  {
    BASE_DIR <- this.path::this.dir()
    setwd(BASE_DIR)
  },
  error = function(e) {
    message("Could not set directory using this.path::this.dir(). Using getwd().")
    BASE_DIR <- getwd()
  }
)

DATA_DIR   <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR   <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR   <- file.path(BASE_DIR, "stan")

dir.create(DATA_DIR,   FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR,   FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP  <- as.integer(Sys.getenv("START_REP", "1"))
message(
  ">>> Pipeline Start (Gamma DGP). Resume settings: Condition=", START_COND,
  ", Replication=", START_REP
)

## -- toggles & constants -------------------------------------------------
RUN_SIM           <- TRUE
RUN_CHECKS        <- TRUE
RUN_FITTING       <- TRUE
RUN_ANALYSIS      <- TRUE
RUN_VISUALIZATION <- FALSE

REPS_PER_CELL <- as.integer(Sys.getenv("REPS_PER_CELL", "200"))
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)

# Global seed is used only for incidental randomness (not the per-cell simulation).
SEED_BASE_SIM <- 2026L
set.seed(SEED_BASE_SIM)

## -- Gamma shape levels (skewness magnitude) ----------------------------
# Comma-separated list, e.g. "2,1,0.5".
# Default keeps k=2 (your current setting) as the least skewed level.
GAMMA_SHAPES_RAW <- Sys.getenv("GAMMA_SHAPES", "2,1,0.5")
GAMMA_SHAPES <- suppressWarnings(as.numeric(strsplit(GAMMA_SHAPES_RAW, ",")[[1]]))
GAMMA_SHAPES <- GAMMA_SHAPES[is.finite(GAMMA_SHAPES) & GAMMA_SHAPES > 0]
GAMMA_SHAPES <- unique(GAMMA_SHAPES)

if (length(GAMMA_SHAPES) < 1) {
  stop("GAMMA_SHAPES must contain at least one positive number (e.g., '2,1,0.5').")
}

# Sort so that the FIRST level is the least-skewed (largest shape).
# This also preserves condition_id ordering in a sensible way.
GAMMA_SHAPES <- sort(GAMMA_SHAPES, decreasing = TRUE)

message(">>> Gamma shape levels (least → most skew): ", paste(GAMMA_SHAPES, collapse = ", "))

## =======================================================================
## 1 · Design grid --------------------------------------------------------
## =======================================================================
var_sets <- list(
  A = matrix(c(0.40, 0.10, 0.10, 0.40), 2, 2, byrow = TRUE),
  B = matrix(c(0.55, 0.10, 0.10, 0.25), 2, 2, byrow = TRUE)
)

assign_distribution_gamma <- function(dir_flag, shape, rate) {
  mirror1 <- substr(dir_flag, 1, 1) == "-"
  mirror2 <- substr(dir_flag, 2, 2) == "-"

  list(
    margin1 = list(type = "gamma", shape = shape, rate = rate, mirror = mirror1),
    margin2 = list(type = "gamma", shape = shape, rate = rate, mirror = mirror2)
  )
}

# Base Study-2 factorial grid (without the Gamma-shape level).
base_grid <- expand.grid(
  direction = c("++", "--", "+-"),
  T         = c(50, 100, 200),
  rho       = c(0.30, 0.50),
  VARset    = names(var_sets),
  stringsAsFactors = FALSE
)

# Add THREE skewness-magnitude levels via shape k.
# NOTE: Ordering is by descending k so the first block corresponds to the
# least-skewed level (by default: k=2), which keeps it aligned with the
# original single-k version.

design_grid <- purrr::map_dfr(GAMMA_SHAPES, function(k) {
  base_grid |>
    mutate(
      shape_k  = k,
      rate_k   = sqrt(k),
      # Use skew_level as the DGP label (Study 1/2 scripts standardize on this name)
      skew_level = paste0("Gamma(k=", k, ")")
    )
}) |>
  mutate(
    condition_id = row_number(),
    n_reps = REPS_PER_CELL
  ) |>
  rowwise() |>
  mutate(
    margin_info = list(assign_distribution_gamma(direction, shape_k, rate_k)),
    phi_matrix  = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid, file = file.path(DATA_DIR, "sim_conditions.rds"))
message(">>> Design grid (Gamma) created (", nrow(design_grid), " conditions).")

## =======================================================================
## 2 · Simulation --------------------------------------------------------
## =======================================================================
if (RUN_SIM) {
  message("\n>>> Running Simulation...")
  source("simulate_data.R", local = TRUE)

  simulate_all_conditions_var1(
    sim_conditions_df = design_grid,
    output_dir        = DATA_DIR,
    start_condition   = START_COND,
    start_rep         = START_REP,
    seed_base         = SEED_BASE_SIM
  )
}

## =======================================================================
## 3 · Quick checks -------------------------------------------------------
## =======================================================================
if (RUN_CHECKS) {
  message("\n>>> Running Simulation Checks...")
  source("check_simulations.R", local = TRUE)
  if (exists("run_post_sim_checks_var1")) {
    run_post_sim_checks_var1(DATA_DIR, CHECKS_DIR)
  }
}

## =======================================================================
## 4 · Stan fitting ------------------------------------------------------
## =======================================================================
if (RUN_FITTING) {
  message("\n>>> Running Stan Fitting...")
  source("fit_models.R", local = TRUE)

  if (exists("fit_var1_copula_models")) {
    fit_var1_copula_models(
      data_dir = DATA_DIR,
      fits_dir = FITS_DIR,
      stan_dir = STAN_DIR,
      results_dir = RESULT_DIR,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15,
      cores_outer = NUM_CORES_OUT,
      start_condition = START_COND,
      start_rep = START_REP
    )
  }
}

## =======================================================================
## 5 · Analysis ----------------------------------------------------------
## =======================================================================
if (RUN_ANALYSIS) {
  message("\n>>> Running Analysis...")
  tryCatch(
    {
      source("analysis_singlelevel.R", local = TRUE)
    },
    error = function(e) {
      message("!!! Analysis failed: ", e$message)
    }
  )
}

## =======================================================================
## 6 · Visualization -----------------------------------------------------
## =======================================================================
if (RUN_VISUALIZATION) {
  message("\n>>> Running Visualization...")
  message("!!! Visualization skipped (toggle RUN_VISUALIZATION to TRUE)")
}

message("\n>>> Pipeline COMPLETE.")
