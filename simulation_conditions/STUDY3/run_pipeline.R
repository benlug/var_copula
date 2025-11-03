###########################################################################
# run_pipeline.R – Study III: Normal margins + Clayton copula
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  if (!requireNamespace("this.path", quietly = TRUE)) {
    if (interactive()) install.packages("this.path")
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

# Study 3 lives in its own folder; no _s3 suffixes
DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR <- file.path(BASE_DIR, "stan")

dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP <- as.integer(Sys.getenv("START_REP", "1"))
message(
  ">>> Pipeline Start (Study III). Resume settings: Condition=", START_COND,
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
set.seed(3033) # dedicated seed for Study 3

## =======================================================================
## 1 · Design grid (Study III)
## =======================================================================
# Only VAR set A (as requested)
var_sets <- list(
  A = matrix(c(
    0.40, 0.10,
    0.10, 0.40
  ), 2, 2, byrow = TRUE)
)

# Normal margins for both series
assign_distribution <- function() {
  list(
    margin1 = list(type = "normal"),
    margin2 = list(type = "normal")
  )
}

# Choose your Clayton θ grid here
THETA_SET <- c(0.5, 1.0, 2.0, 4.0, 8.0)

design_grid <- expand.grid(
  dgp_level = "Normal_Clayton",
  direction = "n/a", # kept for facet compatibility
  T = c(50, 100, 200),
  theta = THETA_SET, # Clayton parameter
  VARset = "A", # ONLY A
  stringsAsFactors = FALSE
) |>
  mutate(
    condition_id = dplyr::row_number(),
    n_reps       = REPS_PER_CELL
  ) |>
  rowwise() |>
  mutate(
    margin_info = list(assign_distribution()),
    phi_matrix  = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid, file = file.path(DATA_DIR, "sim_conditions.rds"))
message(">>> Design grid (Study III) created (", nrow(design_grid), " conditions).")

## =======================================================================
## 2 · Simulation ---------------------------------------------------------
## =======================================================================
if (RUN_SIM) {
  message("\n>>> Running Simulation (S3)...")
  source("simulate_data.R", local = TRUE)
  simulate_all_conditions_var1(
    sim_conditions_df = design_grid,
    output_dir        = DATA_DIR,
    start_condition   = START_COND,
    start_rep         = START_REP
  )
}

## =======================================================================
## 3 · Quick checks -------------------------------------------------------
## =======================================================================
if (RUN_CHECKS) {
  message("\n>>> Running Simulation Checks...")
  source("check_simulations.R", local = TRUE) # reuse your existing file
  if (exists("run_post_sim_checks_var1")) {
    run_post_sim_checks_var1(DATA_DIR, CHECKS_DIR)
  }
}

## =======================================================================
## 4 · Stan fitting -------------------------------------------------------
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
      iter = 3000,
      warmup = 1500,
      adapt_delta = 0.98,
      max_treedepth = 13,
      cores_outer = NUM_CORES_OUT,
      start_condition = START_COND,
      start_rep = START_REP
    )
  }
}

## =======================================================================
## 5 · Analysis -----------------------------------------------------------
## =======================================================================
if (RUN_ANALYSIS) {
  message("\n>>> Running Analysis...")
  tryCatch(
    {
      source("analysis_singlelevel.R", local = TRUE)
    },
    error = function(e) message("!!! Analysis failed: ", e$message)
  )
}

## =======================================================================
## 6 · Visualization ------------------------------------------------------
## =======================================================================
if (RUN_VISUALIZATION) {
  message("\n>>> Visualization not provided (adapt visualize_results.R if desired)")
}

message("\n>>> Pipeline COMPLETE.")
