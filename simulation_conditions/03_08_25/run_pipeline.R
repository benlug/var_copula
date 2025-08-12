###########################################################################
# run_pipeline.R  – Study II: Exponential Tails
###########################################################################

# Load libraries
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

# Define specific directories for Study 2
DATA_DIR <- file.path(BASE_DIR, "data_S2")
CHECKS_DIR <- file.path(BASE_DIR, "checks_S2")
FITS_DIR <- file.path(BASE_DIR, "fits_S2")
RESULT_DIR <- file.path(BASE_DIR, "results_S2")
STAN_DIR <- file.path(BASE_DIR, "stan") # Assuming Stan models are shared

# Create directories
dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP <- as.integer(Sys.getenv("START_REP", "1"))
message(
  ">>> Pipeline Start (Study II). Resume settings: Condition=", START_COND,
  ", Replication=", START_REP
)

## -- toggles & constants -------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FITTING <- TRUE
RUN_ANALYSIS <- TRUE
RUN_VISUALIZATION <- FALSE # Visualization script requires significant updates

REPS_PER_CELL <- 200
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)
set.seed(2026) # New seed for Study 2

## =======================================================================
## 1 · Design grid (Study II) --------------------------------------------
## =======================================================================
var_sets <- list(
  A = matrix(c(0.40, 0.10, 0.10, 0.40), 2, 2, byrow = TRUE),
  B = matrix(c(0.55, 0.10, 0.10, 0.25), 2, 2, byrow = TRUE)
)

# Helper function for Study 2 distribution assignment
assign_distribution_S2 <- function(dir_flag) {
  # Determine mirroring: '+' (Right skew) -> FALSE; '-' (Left skew) -> TRUE
  mirror1 <- substr(dir_flag, 1, 1) == "-"
  mirror2 <- substr(dir_flag, 2, 2) == "-"

  # Standardized Exponential (rate=1)
  list(
    margin1 = list(type = "exponential", rate = 1, mirror = mirror1),
    margin2 = list(type = "exponential", rate = 1, mirror = mirror2)
  )
}

design_grid <- expand.grid(
  # Renaming skew_level to dgp_level for clarity in Study 2
  dgp_level = c("Exponential"),
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
    margin_info = list(assign_distribution_S2(direction)),
    phi_matrix = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid,
  file = file.path(DATA_DIR, "sim_conditions_S2.rds")
)
message(">>> Design grid (Study II) created (", nrow(design_grid), " conditions).")

## =======================================================================
## 2 · Simulation --------------------------------------------------------
## =======================================================================
if (RUN_SIM) {
  message("\n>>> Running Simulation...")
  # We assume simulate_data.R is updated to handle "exponential" type
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
  # We assume fit_models.R is updated for EG model and skew_direction input
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
      # High settings are important due to the boundary constraints of the exponential model
      adapt_delta = 0.995,
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
  # We assume analysis_singlelevel.R is updated for EG parameters
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
  message("!!! Visualization skipped (requires updates to visualize_results.R for Study 2)")
}

message("\n>>> Pipeline COMPLETE.")
