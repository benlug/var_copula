#!/usr/bin/env Rscript
# run_pipeline_ml.R — Study 4 (Multilevel VAR(1), minimal pooling) pipeline
#
# DESIGN (requested in the latest iteration)
#   - One between-unit heterogeneity level (single tau_mu setting)
#   - Two time-series lengths: T = 50 and T = 100
#   - Exponential marginal innovations (EG model correctly specified)
#   - 200 replications per condition (default; override via REPS_PER_CELL)
#
# Parameters held fixed (Study 4 baseline):
#   N_units  = 40
#   burn-in  = 50
#   rho      = 0.50
#   VAR set A: Phi = [[0.40, 0.10], [0.10, 0.40]]
#   direction = "++" (both margins right-skewed; mirroring supported via '-').
#
# IMPORTANT (resume / reproducibility)
#   - Simulation is seeded per (condition_id, rep_id) inside simulate_data_multilevel.R.
#   - Fits are skipped when the corresponding fit_*.rds exists (unless OVERWRITE=TRUE).
#   - If you change the design, consider deleting stale files in ./data and ./fits,
#     or keep them but ensure the fitter filters to the active design conditions.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

## -- base directory ------------------------------------------------------
BASE_DIR <- getwd()
if (requireNamespace("this.path", quietly = TRUE)) {
  tmp <- try(this.path::this.dir(), silent = TRUE)
  if (!inherits(tmp, "try-error") && nzchar(tmp)) BASE_DIR <- tmp
} else if (!is.null(sys.frame(1)$ofile)) {
  # Works when running via source("path/to/run_pipeline_ml.R")
  BASE_DIR <- dirname(normalizePath(sys.frame(1)$ofile))
}
setwd(BASE_DIR)
message("BASE_DIR: ", normalizePath(BASE_DIR, winslash = "/", mustWork = FALSE))

## -- folders -------------------------------------------------------------
DATA_DIR   <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR   <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR   <- file.path(BASE_DIR, "stan")

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHECKS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STAN_DIR, showWarnings = FALSE, recursive = TRUE)

## -- ensure Stan files are where the fitter expects them -----------------
stan_needed <- c("model_EG_ml_min.stan", "model_NG_ml_min.stan")
for (sf in stan_needed) {
  src <- file.path(BASE_DIR, sf)
  dst <- file.path(STAN_DIR, sf)
  if (!file.exists(dst) && file.exists(src)) {
    ok <- file.copy(src, dst, overwrite = TRUE)
    if (isTRUE(ok)) message("Copied ", sf, " -> ", dst)
  }
}
missing_stan <- stan_needed[!file.exists(file.path(STAN_DIR, stan_needed))]
if (length(missing_stan)) {
  stop(
    "Missing Stan model files in: ", STAN_DIR, "\n",
    "Expected: ", paste(missing_stan, collapse = ", "), "\n",
    "Either place the files in 'stan/' or in BASE_DIR so they can be copied."
  )
}

## -- toggles & constants -------------------------------------------------
RUN_SIM    <- as.logical(Sys.getenv("RUN_SIM", "TRUE"))
RUN_CHECKS <- as.logical(Sys.getenv("RUN_CHECKS", "TRUE"))
RUN_FIT    <- as.logical(Sys.getenv("RUN_FIT", "TRUE"))
RUN_ANAL   <- as.logical(Sys.getenv("RUN_ANAL", "TRUE"))

REPS_PER_CELL <- as.integer(Sys.getenv("REPS_PER_CELL", "200"))
START_COND    <- as.integer(Sys.getenv("START_COND", "1"))
START_REP     <- as.integer(Sys.getenv("START_REP", "1"))
OVERWRITE     <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))

CORES_OUTER <- as.integer(Sys.getenv("CORES_OUTER", "0"))
if (is.na(CORES_OUTER) || CORES_OUTER <= 0) {
  CORES_OUTER <- if (requireNamespace("parallel", quietly = TRUE)) {
    max(1L, parallel::detectCores(logical = TRUE) - 1L)
  } else {
    1L
  }
}
message("Outer parallel cores: ", CORES_OUTER)

# Design-level reproducibility (simulation is replication-seeded).
set.seed(2026)

## -- design parameters (shared) -----------------------------------------
N_UNITS <- 40L
BURN_IN <- 50L
rho_val <- 0.50

Phi_A <- matrix(
  c(
    0.40, 0.10,
    0.10, 0.40
  ),
  2, 2, byrow = TRUE
)

direction <- "++"

# Fixed tau_mu (default: Study 4 low-heterogeneity setting).
# Override from the shell, e.g.: TAU_MU="0.40,0.40" Rscript run_pipeline_ml.R
	tau_mu_fixed <- c(0.15, 0.15)
	tau_env <- Sys.getenv("TAU_MU", "")
	if (nzchar(tau_env)) {
	  tau_parsed <- suppressWarnings(as.numeric(strsplit(gsub("[;]", ",", tau_env), "[ ,]+")[[1]]))
	  tau_parsed <- tau_parsed[is.finite(tau_parsed)]
	  if (length(tau_parsed) >= 2) {
	    tau_mu_fixed <- tau_parsed[1:2]
	  } else {
	    warning("Could not parse TAU_MU='", tau_env, "'. Using default (0.15,0.15).")
	  }
	}

assign_distribution_exp <- function(dir_flag) {
  mir1 <- substr(dir_flag, 1, 1) == "-"
  mir2 <- substr(dir_flag, 2, 2) == "-"
  list(
    margin1 = list(type = "exponential", rate = 1, mirror = mir1),
    margin2 = list(type = "exponential", rate = 1, mirror = mir2)
  )
}

T_SER <- c(50L, 100L)

design_ml <- tibble::tibble(
  condition_id = seq_along(T_SER),
  condition    = paste0("C", seq_along(T_SER), "_T", T_SER, "_Exp"),
  level        = "ML_min",
  dgp_level    = "Exponential",
  direction    = direction,
  VARset       = "A",
  T            = T_SER,
  N_units      = N_UNITS,
  burnin       = BURN_IN,
  rho          = rho_val,
  phi_matrix   = rep(list(Phi_A), length(T_SER)),
  tau_mu       = rep(list(tau_mu_fixed), length(T_SER)),
  n_reps       = REPS_PER_CELL,
  margin_info  = rep(list(assign_distribution_exp(direction)), length(T_SER))
)

saveRDS(design_ml, file.path(DATA_DIR, "sim_conditions_ml.rds"))
message("\n>>> ML design saved: ", file.path(DATA_DIR, "sim_conditions_ml.rds"))
message("    Conditions: ", nrow(design_ml), " (", paste(design_ml$condition, collapse = ", "), ")")
message("    Replications per condition: ", REPS_PER_CELL)
message(sprintf("    tau_mu fixed at (%.3f, %.3f)", tau_mu_fixed[1], tau_mu_fixed[2]))

## -- simulation ----------------------------------------------------------
if (RUN_SIM) {
  message("\n>>> Running ML simulation …")
  source("simulate_data_multilevel.R", local = TRUE)
  simulate_all_conditions_var1_ml(
    sim_conditions_df = design_ml,
    output_dir        = DATA_DIR,
    start_condition   = START_COND,
    start_rep         = START_REP,
    overwrite         = OVERWRITE
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

  if (!file.exists(file.path(BASE_DIR, "fit_models_ml.R"))) {
    stop(
      "Missing 'fit_models_ml.R' in BASE_DIR.\n",
      "Expected at: ", file.path(BASE_DIR, "fit_models_ml.R")
    )
  }

  source("fit_models_ml.R", local = TRUE)

  fit_var1_copula_models_ml(
    data_dir = DATA_DIR,
    fits_dir = FITS_DIR,
    stan_dir = STAN_DIR,
    results_dir = RESULT_DIR,
    chains = as.integer(Sys.getenv("CHAINS", "4")),
    iter = as.integer(Sys.getenv("ITER", "3000")),
    warmup = as.integer(Sys.getenv("WARMUP", "1500")),
    adapt_delta = as.numeric(Sys.getenv("ADAPT_DELTA", "0.90")),
    max_treedepth = as.integer(Sys.getenv("MAX_TD", "10")),
    cores_outer = CORES_OUTER,
    start_condition = START_COND,
    start_rep = START_REP,
    overwrite = OVERWRITE
  )
}

## -- analysis ------------------------------------------------------------
if (RUN_ANAL) {
  message("\n>>> Analyzing ML results …")
  source("analysis_multilevel.R", local = TRUE)
}

message("\n>>> ML pipeline COMPLETE.")
