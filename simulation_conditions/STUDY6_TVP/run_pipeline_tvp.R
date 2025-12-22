###########################################################################
# run_pipeline_tvp.R  – Study 6: Time-Varying Parameter Copula-VAR
###########################################################################
#
# Research Questions:
#   1. Detection power: Under what conditions can we detect time-varying rho?
#   2. Bias from ignoring TVP: How biased are constant-rho estimates when rho varies?
#   3. False positive rate: When rho is constant, how often do TVP models suggest variation?
#   4. Computational feasibility: Can models run on T=50-200 with reasonable diagnostics?
#
###########################################################################

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  if (!requireNamespace("this.path", quietly = TRUE)) {
    if (interactive()) install.packages("this.path")
    library(this.path)
  }
})

## -- folders -----------------------------------------------------------------
tryCatch({
  BASE_DIR <- this.path::this.dir()
  setwd(BASE_DIR)
}, error = function(e) {
  message("Could not set directory using this.path::this.dir(). Using getwd().")
  BASE_DIR <- getwd()
})

DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR <- file.path(BASE_DIR, "stan")

dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- resume flags ------------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP <- as.integer(Sys.getenv("START_REP", "1"))
message(">>> Pipeline Start (Study 6 - TVP). Resume: Condition=", START_COND,
        ", Replication=", START_REP)

## -- toggles & constants -----------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- FALSE
RUN_FITTING <- TRUE
RUN_ANALYSIS <- TRUE

REPS_PER_CELL <- 200
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)
set.seed(2028)  # Study 6 seed

## ==========================================================================
## 1. Design grid (Study 6: TVP-Copula-VAR)
## ==========================================================================

var_sets <- list(
  A = matrix(c(0.40, 0.10, 0.10, 0.40), 2, 2, byrow = TRUE)
)

# Define TVP patterns and their parameters
define_tvp_patterns <- function() {
  list(
    # Pattern 1: Constant (null case)
    constant = list(
      pattern = "constant",
      params = list(rho_base = 0.5)
    ),
    
    # Pattern 2: Linear drift (therapy effect)
    linear_down = list(
      pattern = "linear",
      params = list(rho_start = 0.7, rho_end = 0.2)
    ),
    
    # Pattern 3: Step change (episode onset)
    step_up = list(
      pattern = "step",
      params = list(rho_before = 0.3, rho_after = 0.7, tau = 0.5)
    ),
    
    # Pattern 4: Random walk (general TVP)
    random_walk_moderate = list(
      pattern = "random_walk",
      params = list(rho_start = 0.5, sigma_z = 0.05)
    )
  )
}

# Define margin configurations
define_margins <- function(margin_type, direction) {
  mirror1 <- substr(direction, 1, 1) == "-"
  mirror2 <- substr(direction, 2, 2) == "-"
  
  if (margin_type == "normal") {
    list(
      margin1 = list(type = "normal", mirror = FALSE),
      margin2 = list(type = "normal", mirror = FALSE)
    )
  } else {
    list(
      margin1 = list(type = "exponential", rate = 1, mirror = mirror1),
      margin2 = list(type = "exponential", rate = 1, mirror = mirror2)
    )
  }
}

# Build design grid
tvp_patterns <- define_tvp_patterns()

design_grid <- expand.grid(
  T = c(100, 200),
  tvp_pattern_name = names(tvp_patterns),
  margin_type = c("normal", "exponential"),
  direction = c("++"),
  VARset = names(var_sets),
  stringsAsFactors = FALSE
) |>
  mutate(
    condition_id = row_number(),
    n_reps = REPS_PER_CELL
  ) |>
  rowwise() |>
  mutate(
    tvp_pattern = tvp_patterns[[tvp_pattern_name]]$pattern,
    tvp_params = list(tvp_patterns[[tvp_pattern_name]]$params),
    margin_info = list(define_margins(margin_type, direction)),
    phi_matrix = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid, file = file.path(DATA_DIR, "sim_conditions.rds"))
message(">>> Design grid created (", nrow(design_grid), " conditions x ", 
        REPS_PER_CELL, " reps = ", nrow(design_grid) * REPS_PER_CELL, " fits per model).")

# Print design summary
message("\n>>> Design Summary:")
design_grid |>
  group_by(T, tvp_pattern, margin_type) |>
  summarise(n_conditions = n(), .groups = "drop") |>
  print()

## ==========================================================================
## 2. Simulation
## ==========================================================================
if (RUN_SIM) {
  message("\n>>> Running TVP Simulation...")
  source("simulate_data_tvp.R", local = TRUE)
  
  simulate_all_conditions_tvp(
    sim_conditions_df = design_grid,
    output_dir = DATA_DIR,
    start_condition = START_COND,
    start_rep = START_REP
  )
}

## ==========================================================================
## 3. Checks (optional)
## ==========================================================================
if (RUN_CHECKS) {
  message("\n>>> Running Simulation Checks...")
  # Add check functions here if needed
}

## ==========================================================================
## 4. Stan fitting
## ==========================================================================
if (RUN_FITTING) {
  message("\n>>> Running Stan Fitting...")
  source("fit_models_tvp.R", local = TRUE)
  
  if (exists("fit_tvp_copula_models")) {
    fit_tvp_copula_models(
      data_dir = DATA_DIR,
      fits_dir = FITS_DIR,
      stan_dir = STAN_DIR,
      results_dir = RESULT_DIR,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.9,
      max_treedepth = 12,
      cores_outer = NUM_CORES_OUT,
      start_condition = START_COND,
      start_rep = START_REP,
      models_to_fit = c("TVP_NG", "Const_NG", "TVP_EG", "Const_EG")
    )
  }
}

## ==========================================================================
## 5. Analysis
## ==========================================================================
if (RUN_ANALYSIS) {
  message("\n>>> Running Analysis...")
  tryCatch({
    source("analysis_tvp.R", local = TRUE)
  }, error = function(e) {
    message("!!! Analysis failed: ", e$message)
  })
}

message("\n>>> Pipeline COMPLETE.")
