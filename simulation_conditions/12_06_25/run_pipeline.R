###########################################################################
# run_pipeline.R  – Single‑Level VAR(1) Copula simulation (72 cells)
###########################################################################

library(dplyr)
library(purrr)
library(tidyr)
library(this.path)

## -- folders -------------------------------------------------------------
BASE_DIR <- this.dir()
setwd(BASE_DIR)
DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULT_DIR <- file.path(BASE_DIR, "results")
dir.create(DATA_DIR, FALSE, TRUE)
dir.create(CHECKS_DIR, FALSE, TRUE)
dir.create(FITS_DIR, FALSE, TRUE)
dir.create(RESULT_DIR, FALSE, TRUE)

## -- resume flags --------------------------------------------------------
START_COND <- as.integer(Sys.getenv("START_COND", "1")) # 1‑based index
START_REP <- as.integer(Sys.getenv("START_REP", "1")) # within‑cell
message(
  ">>> Resume settings: start at condition ", START_COND,
  ", replication ", START_REP
)

## -- toggles & constants -------------------------------------------------
RUN_SIM <- TRUE
RUN_CHECKS <- TRUE
RUN_FITTING <- TRUE
REPS_PER_CELL <- 100
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
  T = c(50, 100),
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

## =======================================================================
## 2 · Simulation  (skips finished cells, supports resume) ---------------
## =======================================================================
if (RUN_SIM) {
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
  source("check_simulations.R", local = TRUE)
  run_post_sim_checks_var1(DATA_DIR, CHECKS_DIR)
}

## =======================================================================
## 4 · Stan fitting (resume‑aware) ---------------------------------------
## =======================================================================
if (RUN_FITTING) {
  source("fit_models.R", local = TRUE)
  fit_var1_copula_models(
    data_dir = DATA_DIR,
    fits_dir = FITS_DIR,
    stan_dir = file.path(BASE_DIR, "stan"),
    results_dir = RESULT_DIR,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    cores_outer = NUM_CORES_OUT,
    start_condition = START_COND,
    start_rep = START_REP
  )
}

message("\n>>> Pipeline COMPLETE.")
