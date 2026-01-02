###########################################################################
# run_pipeline.R  - Single-Level VAR(1) Copula Simulation Pipeline
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  if (!requireNamespace("this.path", quietly = TRUE)) {
    if (interactive()) {
      install.packages("this.path")
    }
  }
  if (!requireNamespace("this.path", quietly = TRUE)) {
    stop("Package 'this.path' is required to locate the Study 1 directory. Please install it.")
  }
  library(this.path)
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
    setwd(BASE_DIR)
  }
)
message("Working directory: ", BASE_DIR)

DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
STAN_DIR <- file.path(BASE_DIR, "stan")

dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)
dir.create(STAN_DIR, FALSE, TRUE)

## -- Verify Stan files exist ----------------------------------------------
stan_files <- c(
  file.path(STAN_DIR, "model_SNG_sl.stan"),
  file.path(STAN_DIR, "model_NG_sl.stan")
)
missing_stan <- stan_files[!file.exists(stan_files)]
if (length(missing_stan) > 0) {
  stop(
    "Missing Stan files:\n  ", paste(missing_stan, collapse = "\n  "),
    "\nPlease place Stan model files in: ", STAN_DIR
  )
}

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1"))
START_REP <- as.integer(Sys.getenv("START_REP", "1"))
message(">>> Pipeline Start. Resume: Condition=", START_COND, ", Rep=", START_REP)

## -- toggles -------------------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FITTING <- TRUE
RUN_ANALYSIS <- TRUE
RUN_VISUALIZATION <- TRUE

REPS_PER_CELL <- 200
BURNIN <- 100
NUM_CORES_OUT <- max(1, parallel::detectCores() - 10)
set.seed(2025)

## =======================================================================
## 1. Design Grid
## =======================================================================
var_sets <- list(
  A = matrix(c(
    0.40, 0.10,
    0.10, 0.40
  ), 2, 2, byrow = TRUE),
  B = matrix(c(
    0.55, 0.10,
    0.10, 0.25
  ), 2, 2, byrow = TRUE)
)

assign_skew <- function(level, dir_flag) {
  sgn <- function(x) ifelse(x == "+", 1, -1)
  s1 <- sgn(substr(dir_flag, 1, 1))
  s2 <- sgn(substr(dir_flag, 2, 2))

  switch(level,
    moderateSN = list(
      margin1 = list(type = "skewnormal", alpha = s1 * 4),
      margin2 = list(type = "skewnormal", alpha = s2 * 4)
    ),
    strongSN = list(
      margin1 = list(type = "skewnormal", alpha = s1 * 9),
      margin2 = list(type = "skewnormal", alpha = s2 * 9)
    ),
    extremeCHI = list(
      margin1 = list(type = "chisq", df = 1, mirror = (s1 < 0)),
      margin2 = list(type = "chisq", df = 1, mirror = (s2 < 0))
    ),
    stop("Unknown skew level: ", level)
  )
}

design_grid <- expand.grid(
  skew_level = c("moderateSN", "strongSN", "extremeCHI"),
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
    margin_info = list(assign_skew(skew_level, direction)),
    phi_matrix  = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid, file.path(DATA_DIR, "sim_conditions_singlelevel.rds"))
message(">>> Design grid: ", nrow(design_grid), " conditions x ", REPS_PER_CELL, " reps")

## =======================================================================
## 2. Simulation
## =======================================================================
if (RUN_SIM) {
  message("\n>>> Running Simulation (burn-in = ", BURNIN, ")...")
  source("simulate_data.R", local = TRUE)

  simulate_all_conditions_var1(
    sim_conditions_df = design_grid,
    output_dir        = DATA_DIR,
    burnin            = BURNIN,
    start_condition   = START_COND,
    start_rep         = START_REP,
    verify_copula     = FALSE
  )
}

## =======================================================================
## 3. Simulation Checks
## =======================================================================
if (RUN_CHECKS) {
  message("\n>>> Running Simulation Checks...")
  source("check_simulations.R", local = TRUE)
  if (exists("run_post_sim_checks_var1")) {
    run_post_sim_checks_var1(DATA_DIR, CHECKS_DIR)
  }
}

## =======================================================================
## 4. Stan Fitting
## =======================================================================
if (RUN_FITTING) {
  message("\n>>> Running Stan Fitting...")
  source("fit_models.R", local = TRUE)

  if (exists("fit_var1_copula_models")) {
    fit_var1_copula_models(
      data_dir        = DATA_DIR,
      fits_dir        = FITS_DIR,
      stan_dir        = STAN_DIR,
      results_dir     = RESULT_DIR,
      chains          = 4,
      iter            = 4000,
      warmup          = 2000,
      adapt_delta     = 0.95,
      max_treedepth   = 15,
      cores_outer     = NUM_CORES_OUT,
      start_condition = START_COND,
      start_rep       = START_REP
    )
  }
}

## =======================================================================
## 5. Analysis
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
## 6. Visualization
## =======================================================================
if (RUN_VISUALIZATION) {
  message("\n>>> Running Visualization...")
  tryCatch(
    {
      source("visualize_results.R", local = TRUE)
    },
    error = function(e) {
      message("!!! Visualization failed: ", e$message)
    }
  )
}

message("\n>>> Pipeline COMPLETE.")
