###########################################################################
# run_pipeline.R  – Single‑Level VAR(1) Copula simulation
###########################################################################

# Load libraries; prefer using suppressPackageStartupMessages for cleaner logs
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  # Ensure this.path is available for robust directory setting
  if (!requireNamespace("this.path", quietly = TRUE)) {
    # If running this script interactively, attempt to install if missing
    if (interactive()) {
      print("Installing 'this.path' package for robust directory handling.")
      install.packages("this.path")
    }
    library(this.path)
  }
})


## -- folders -------------------------------------------------------------
# Robustly set the base directory to the script's location
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


DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
# Added Stan directory definition
STAN_DIR <- file.path(BASE_DIR, "stan")

# Create directories recursively and silently
dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)
# dir.create(STAN_DIR, FALSE, TRUE) # Assuming stan dir already exists

## -- resume flags --------------------------------------------------------
# Allow resuming via environment variables (useful for job arrays)
START_COND <- as.integer(Sys.getenv("START_COND", "1")) # 1‑based index
START_REP <- as.integer(Sys.getenv("START_REP", "1")) # within‑cell
message(
  ">>> Pipeline Start. Resume settings: Condition=", START_COND,
  ", Replication=", START_REP
)

## -- toggles & constants -------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FITTING <- TRUE
# Added toggles for analysis and visualization steps
RUN_ANALYSIS <- TRUE
RUN_VISUALIZATION <- TRUE

REPS_PER_CELL <- 200
# Use available cores wisely, leaving one free for system responsiveness
NUM_CORES_OUT <- max(1, parallel::detectCores() - 1)
set.seed(2025)

## =======================================================================
## 1 · Design grid -------------------------------------------------------
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
    # Chi-squared df=1 is highly skewed. Mirroring provides left skew.
    extremeCHI = list(
      margin1 = list(type = "chisq", df = 1, mirror = (s1 < 0)),
      margin2 = list(type = "chisq", df = 1, mirror = (s2 < 0))
    ),
    stop("unknown skew level")
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
    phi_matrix = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid,
  file = file.path(DATA_DIR, "sim_conditions_singlelevel.rds")
)
message(">>> Design grid created (", nrow(design_grid), " conditions).")

## =======================================================================
## 2 · Simulation  (skips finished cells, supports resume) ---------------
## =======================================================================
if (RUN_SIM) {
  message("\n>>> Running Simulation...")
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
  # Check if the function exists before calling it (robustness)
  if (exists("run_post_sim_checks_var1")) {
    run_post_sim_checks_var1(DATA_DIR, CHECKS_DIR)
  }
}

## =======================================================================
## 4 · Stan fitting (resume‑aware) ---------------------------------------
## =======================================================================
if (RUN_FITTING) {
  message("\n>>> Running Stan Fitting...")
  source("fit_models.R", local = TRUE)

  # Stan settings: Increased adapt_delta and max_treedepth for robust fitting
  # These settings are generally recommended for complex models.
  if (exists("fit_var1_copula_models")) {
    fit_var1_copula_models(
      data_dir = DATA_DIR,
      fits_dir = FITS_DIR,
      stan_dir = STAN_DIR,
      results_dir = RESULT_DIR,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      adapt_delta = 0.99, # High adapt_delta to minimize divergences
      max_treedepth = 15, # Increased depth to handle complex geometry
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
  # Sourcing the analysis script; execution happens within the script.
  # Use tryCatch to prevent a failure here from stopping visualization
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
