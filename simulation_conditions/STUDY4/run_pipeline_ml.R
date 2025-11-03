#!/usr/bin/env Rscript
# run_pipeline_ml.R — Multilevel (minimal pooling) pilot: N=50, T=100, burn-in=30
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
    message("this.path::this.dir() failed; using getwd().")
    BASE_DIR <<- getwd()
  }
)
DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR <- file.path(BASE_DIR, "stan")
dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- toggles & constants -------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FIT <- TRUE
RUN_ANAL <- TRUE

REPS_PER_CELL <- as.integer(Sys.getenv("REPS_PER_CELL", "50")) # pilot size
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)
set.seed(2026)

## -- pilot design: 2 conditions -----------------------------------------
# fixed pieces
N_UNITS <- 50L
T_SER <- 100L
BURN_IN <- 30L
rho_val <- 0.50
Phi_A <- matrix(c(0.40, 0.10, 0.10, 0.40), 2, 2, byrow = TRUE) # always stable here
direction <- "++"

assign_distribution <- function(dir_flag) {
  mir1 <- substr(dir_flag, 1, 1) == "-"
  mir2 <- substr(dir_flag, 2, 2) == "-"
  list(
    margin1 = list(type = "exponential", rate = 1, mirror = mir1),
    margin2 = list(type = "exponential", rate = 1, mirror = mir2)
  )
}

design_ml <- tibble::tibble(
  condition_id = 1:2,
  condition = c("C1_lowHet", "C2_highHet"),
  level = "ML_min",
  dgp_level = "Exponential",
  direction = direction,
  VARset = "A",
  T = T_SER,
  N_units = N_UNITS,
  burnin = BURN_IN,
  rho = rho_val,
  phi_matrix = list(Phi_A, Phi_A),
  tau_mu = list(c(0.15, 0.15), c(0.40, 0.40)), # only difference
  n_reps = REPS_PER_CELL,
  margin_info = list(
    assign_distribution(direction),
    assign_distribution(direction)
  )
)

saveRDS(design_ml, file.path(DATA_DIR, "sim_conditions_ml.rds"))
message(">>> ML design created with ", nrow(design_ml), " conditions.")

## -- simulation ----------------------------------------------------------
if (RUN_SIM) {
  message("\n>>> Running ML simulation …")
  source("simulate_data_multilevel.R", local = TRUE)
  simulate_all_conditions_var1_ml(
    sim_conditions_df = design_ml,
    output_dir        = DATA_DIR,
    start_condition   = as.integer(Sys.getenv("START_COND", "1")),
    start_rep         = as.integer(Sys.getenv("START_REP", "1")),
    overwrite         = as.logical(Sys.getenv("OVERWRITE", "FALSE"))
  )
}

## -- checks --------------------------------------------------------------
if (RUN_CHECKS) {
  message("\n>>> Running quick ML checks …")
  source("check_simulations_ml.R", local = TRUE)
  if (exists("run_post_sim_checks_var1_ml")) {
    run_post_sim_checks_var1_ml(DATA_DIR, CHECKS_DIR, n_units_to_plot = 3)
  }
}

## -- fitting -------------------------------------------------------------
if (RUN_FIT) {
  message("\n>>> Fitting ML models (EG_MLmin & NG_MLmin) …")
  source("fit_models_ml.R", local = TRUE)
  fit_var1_copula_models_ml(
    data_dir = DATA_DIR,
    fits_dir = FITS_DIR,
    stan_dir = STAN_DIR,
    results_dir = RESULT_DIR,
    chains = 4, iter = 4000, warmup = 2000,
    adapt_delta = 0.997, max_treedepth = 15,
    cores_outer = NUM_CORES_OUT,
    start_condition = as.integer(Sys.getenv("START_COND", "1")),
    start_rep = as.integer(Sys.getenv("START_REP", "1"))
  )
}

## -- analysis ------------------------------------------------------------
if (RUN_ANAL) {
  message("\n>>> Analyzing ML results …")
  source("analysis_multilevel.R", local = TRUE)
}
message("\n>>> ML pipeline COMPLETE.")
