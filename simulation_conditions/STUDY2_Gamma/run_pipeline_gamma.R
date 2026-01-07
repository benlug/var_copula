#!/usr/bin/env Rscript
###########################################################################
# run_pipeline_gamma.R  – Gamma Innovations (Study 3)
#
# Same factorial design as Study 2, except the DGP innovations use
# standardized Gamma margins instead of standardized Exponential margins.
#
# Notes
#   * We use Gamma(shape = k, rate = sqrt(k)) so that Var(X)=1.
#     Then Z = X - E[X] = X - sqrt(k) has mean 0 and variance 1.
#   * Mirroring ("-") multiplies Z by -1 to create left-skewness.
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

# Define specific directories
DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR <- file.path(BASE_DIR, "stan")

# Create directories
dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP <- as.integer(Sys.getenv("START_REP", "1"))
message(
  ">>> Pipeline Start (Gamma DGP). Resume settings: Condition=", START_COND,
  ", Replication=", START_REP
)

## -- toggles & constants -------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FITTING <- TRUE
RUN_ANALYSIS <- TRUE
RUN_VISUALIZATION <- FALSE

REPS_PER_CELL <- 200
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)

# Fixed Gamma shape (k). Keep it constant to preserve the Study 2 factorial design.
GAMMA_SHAPE <- suppressWarnings(as.numeric(Sys.getenv("GAMMA_SHAPE", "2")))
if (!is.finite(GAMMA_SHAPE) || GAMMA_SHAPE <= 0) {
  stop("GAMMA_SHAPE must be a positive number")
}
GAMMA_RATE <- sqrt(GAMMA_SHAPE)

# Global seed is only used for non-simulation randomness.
SEED_BASE_SIM <- 2026L
set.seed(SEED_BASE_SIM)

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

design_grid <- expand.grid(
  dgp_level = c("Gamma"),
  direction = c("++", "--", "+-"),
  T = c(50, 100, 200),
  rho = c(0.30, 0.50),
  VARset = names(var_sets),
  stringsAsFactors = FALSE
) |>
  mutate(
    condition_id = row_number(),
    n_reps = REPS_PER_CELL
  ) |>
  rowwise() |>
  mutate(
    margin_info = list(assign_distribution_gamma(direction, GAMMA_SHAPE, GAMMA_RATE)),
    phi_matrix = list(var_sets[[VARset]])
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
    output_dir = DATA_DIR,
    start_condition = START_COND,
    start_rep = START_REP,
    seed_base = SEED_BASE_SIM
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
