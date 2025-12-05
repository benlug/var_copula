# run_pipeline_SEM.R — Study 5 : SEM Skewness (A: indicator, B: latent)
# Pipe: simulate → (optional) checks → fit → analyze
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
    message("Could not set directory with this.path; using getwd().")
    BASE_DIR <- getwd()
  }
)

DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks_sem")
FITS_DIR <- file.path(BASE_DIR, "fits_sem")
RESULT_DIR <- file.path(BASE_DIR, "results_sem")
STAN_DIR <- file.path(BASE_DIR, "stan")

dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)
dir.create(STAN_DIR, FALSE, TRUE)

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP <- as.integer(Sys.getenv("START_REP", "1"))
message(">>> SEM Pipeline Start. Condition=", START_COND, "  Replication=", START_REP)

## -- toggles & constants -------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FITTING <- TRUE
RUN_ANALYSIS <- TRUE
RUN_VISUALIZATION <- FALSE

REPS_PER_CELL <- 90
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)
set.seed(2042)

## =======================================================================
## Design grid ------------------------------------------------------------
## =======================================================================
B_mat <- matrix(c(
  0.55, 0.10,
  0.10, 0.25
), 2, 2, byrow = TRUE)

assign_skew <- function(dir_flag) {
  s1 <- if (substr(dir_flag, 1, 1) == "+") +1L else -1L
  s2 <- if (substr(dir_flag, 2, 2) == "+") +1L else -1L
  c(s1, s2)
}

design_grid <- expand.grid(
  sem_study = c("A_indicator", "B_latent"),
  direction = c("++", "+-"),
  T = c(100, 200),
  rho = c(0.30),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    condition_id = dplyr::row_number(),
    n_reps       = REPS_PER_CELL,
    phi_matrix   = list(B_mat)[rep(1, dplyr::n())],
    skew_dir     = purrr::map(direction, assign_skew)
  )

saveRDS(design_grid, file = file.path(DATA_DIR, "sim_conditions_sem.rds"))
message(">>> SEM design created (", nrow(design_grid), " conditions).")

## =======================================================================
## 1 · Simulation ---------------------------------------------------------
## =======================================================================
if (RUN_SIM) {
  message("\n>>> Running SEM Simulation …")
  source("simulate_data_SEM.R", local = TRUE)
  simulate_all_conditions_SEM(
    sim_conditions_df = design_grid,
    output_dir        = DATA_DIR,
    start_condition   = START_COND,
    start_rep         = START_REP
  )
}

## =======================================================================
## 2 · Checks (DGP visualizations) ----------------------------------------
## =======================================================================
if (RUN_CHECKS) {
  message("\n>>> Running Checks …")
  source("check_simulations_SEM.R", local = TRUE)
  run_post_sim_checks_sem(DATA_DIR, CHECKS_DIR)
}

## =======================================================================
## 3 · Stan fitting -------------------------------------------------------
## =======================================================================
if (RUN_FITTING) {
  message("\n>>> Running Stan Fitting (SEM) …")
  source("fit_models_SEM.R", local = TRUE)
  fit_sem_models(
    data_dir        = DATA_DIR,
    fits_dir        = FITS_DIR,
    stan_dir        = STAN_DIR,
    results_dir     = RESULT_DIR,
    chains          = 4,
    iter            = 4000,
    warmup          = 2000,
    adapt_delta     = 0.95, # per your request
    max_treedepth   = 12, # per your request
    cores_outer     = NUM_CORES_OUT,
    start_condition = START_COND,
    start_rep       = START_REP
  )
}

## =======================================================================
## 4 · Analysis -----------------------------------------------------------
## =======================================================================
if (RUN_ANALYSIS) {
  message("\n>>> Running SEM Analysis …")
  source("analysis_sem.R", local = TRUE)
}

if (RUN_VISUALIZATION) {
  message("\n>>> Visualization skipped (handled in analysis doc).")
}

message("\n>>> SEM Pipeline COMPLETE.")
