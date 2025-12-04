###########################################################################
# run_pipeline.R  –  Master driver for the multilevel VAR(1) copula project
###########################################################################

library(dplyr)
library(this.path)

# ----------------------------- USER SETTINGS -----------------------------
## Factorial design -------------------------------------------------------
# SKEW_LEVELS <- c("moderateSN", "strongSN", "extremeCHI")
# DIRECTIONS <- c("++", "--", "+-", "-+")
# TIME_POINTS <- c(50, 100)
# RHO_LEVELS <- c(0.30, 0.50)

SKEW_LEVELS <- c("moderateSN")
DIRECTIONS <- c("++")
TIME_POINTS <- c(100)
RHO_LEVELS <- c(0.30)

## Population Φ matrices --------------------------------------------------
VAR_SETS <- list(
  A = matrix(c(0.40, 0.10, 0.10, 0.40), 2, 2, byrow = TRUE)
  # B = matrix(c(0.40, -0.10, -0.10, 0.40), 2, 2, byrow = TRUE),
  # C = matrix(c(0.55, 0.10, 0.10, 0.25), 2, 2, byrow = TRUE),
  # D = matrix(c(0.25, 0.15, 0.05, 0.55), 2, 2, byrow = TRUE)
)

## Simulation design ------------------------------------------------------
REPS_PER_CELL <- 2 # independent repetitions per design cell
N_SUBJECTS <- 100
SD_AR <- 0.10
SD_CL <- 0.05
BURNIN <- 30

## Parallelisation --------------------------------------------------------
CORES_SIM <- parallel::detectCores() - 2 # ← number of cores for data simulation
CORES_FIT <- parallel::detectCores() - 14 # ← number of cores for dataset‑level Stan fits

## Stan settings ----------------------------------------------------------
CHAINS <- 4
PARALLEL_CHAINS <- TRUE # TRUE ⇒ options(mc.cores = CHAINS)
STAN_ITER <- 4000
STAN_WARMUP <- 2000
STAN_ADAPT_DELTA <- 0.90
STAN_MAX_TREEDPTH <- 12

## Which stages to run ----------------------------------------------------
DO_SIMULATION <- TRUE
DO_CHECKS <- TRUE
DO_FITTING <- TRUE
DO_ANALYSIS <- TRUE
# ------------------------------------------------------------------------

# --------------------------- folder layout ------------------------------
BASE_DIR <- this.dir()
DATA_DIR <- file.path(BASE_DIR, "data")
CHECK_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")
dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECK_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RES_DIR, FALSE, TRUE)

# --------------------------- source helpers -----------------------------
source(file.path(BASE_DIR, "simulate_data.R"))
source(file.path(BASE_DIR, "check_simulations.R"))

# --------------------------- design grid --------------------------------
design_grid <- build_multilevel_design_grid(
  skew_levels = SKEW_LEVELS,
  directions  = DIRECTIONS,
  T_vals      = TIME_POINTS,
  rhos        = RHO_LEVELS,
  var_sets    = VAR_SETS,
  reps        = REPS_PER_CELL
)

# --------------------------- 1. simulation ------------------------------
if (DO_SIMULATION) {
  simulate_all_conditions_ml(
    sim_conditions_df = design_grid,
    output_dir        = DATA_DIR,
    N_per_subject     = N_SUBJECTS,
    burnin            = BURNIN,
    sd_AR             = SD_AR,
    sd_CL             = SD_CL,
    cores             = CORES_SIM # <<< NEW
  )
}

# --------------------------- 2. quick PDFs ------------------------------
if (DO_CHECKS) {
  run_post_sim_checks_ml(
    data_dir           = DATA_DIR,
    checks_dir         = CHECK_DIR,
    n_subjects_to_plot = 3,
    n_reps_to_plot     = 1
  )
}

# --------------------------- 3. Stan fits -------------------------------
if (DO_FITTING) {
  source(file.path(BASE_DIR, "fit_models.R"))
  fit_models_ml(
    data_dir         = DATA_DIR,
    fits_dir         = FITS_DIR,
    stan_dir         = file.path(BASE_DIR, "stan"),
    chains           = CHAINS,
    iter             = STAN_ITER,
    warmup           = STAN_WARMUP,
    adapt_delta      = STAN_ADAPT_DELTA,
    max_treedepth    = STAN_MAX_TREEDPTH,
    parallel_chains  = PARALLEL_CHAINS,
    cores_outer      = CORES_FIT # <<< NEW
  )
}

# --------------------------- 4. analysis --------------------------------
if (DO_ANALYSIS) {
  source(file.path(BASE_DIR, "analyze_results.R"))
  analyze_results_ml(DATA_DIR, FITS_DIR, RES_DIR)
}

message("\n>>> Pipeline COMPLETE.")
