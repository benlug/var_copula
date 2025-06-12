###########################################################################
# run_pipeline.R  –  Single‑Level VAR(1) Copula Simulation (192 cells)
# -------------------------------------------------------------------------
# • 3 skew magnitudes  × 4 direction patterns  × 2 T  × 2 ρ  × 4 VAR sets
# • 100 replications per design cell
# • Fits two models:  SkewNormal+G‑Copula  (correct)   vs.   Normal+G‑Copula
# -------------------------------------------------------------------------
# Requires:
#   simulate_data.R
#   check_simulations.R
#   fit_models.R
#   Stan files: model_SNG_sl.stan   model_NG_sl.stan
###########################################################################

# --- Setup ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(this.path)

BASE_DIR <- this.dir()
setwd(BASE_DIR)

DATA_DIR <- file.path(BASE_DIR, "data")
CHECKS_DIR <- file.path(BASE_DIR, "checks")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULTS_DIR <- file.path(BASE_DIR, "results")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHECKS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Control flags --------------------------------------------------------
RUN_SINGLE_SIM <- TRUE # main 192‑cell simulation
RUN_SINGLE_CHECK <- TRUE
RUN_SINGLE_FIT <- TRUE
RUN_EXP_TEST <- FALSE # will be in a separate script
RUN_MULTILEVEL <- FALSE # phase‑2 scripts later

NUM_CORES <- parallel::detectCores() - 2
REPS_PER_CELL <- NUM_CORES

# -------------------------------------------------------------------------
# 1.  DESIGN GRID (192 conditions) ----------------------------------------
# -------------------------------------------------------------------------
set.seed(2025)

# (a) VAR coefficient matrices
var_sets <- list(
  A = matrix(c(
    0.40, 0.10,
    0.10, 0.40
  ), nrow = 2, byrow = TRUE),
  B = matrix(c(
    0.40, -0.10,
    -0.10, 0.40
  ), nrow = 2, byrow = TRUE),
  C = matrix(c(
    0.55, 0.10,
    0.10, 0.25
  ), nrow = 2, byrow = TRUE),
  D = matrix(c(
    0.25, 0.15,
    0.05, 0.55
  ), nrow = 2, byrow = TRUE)
)

# (b) helper to assign marginal shape info
assign_skew <- function(level, dir_flag) {
  sign1 <- ifelse(substr(dir_flag, 1, 1) == "+", 1, -1)
  sign2 <- ifelse(substr(dir_flag, 2, 2) == "+", 1, -1)
  if (level == "moderateSN") {
    list(
      margin1 = list(type = "skewnormal", alpha = sign1 * 4),
      margin2 = list(type = "skewnormal", alpha = sign2 * 4)
    )
  } else if (level == "strongSN") {
    list(
      margin1 = list(type = "skewnormal", alpha = sign1 * 9),
      margin2 = list(type = "skewnormal", alpha = sign2 * 9)
    )
  } else { # extremeCHI
    list(
      margin1 = list(type = "chisq", df = 1, mirror = (sign1 < 0)),
      margin2 = list(type = "chisq", df = 1, mirror = (sign2 < 0))
    )
  }
}

design_grid <- expand.grid(
  skew_level = c("moderateSN", "strongSN", "extremeCHI"),
  direction = c("++", "--", "+-", "-+"),
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
    phi_matrix  = list(var_sets[[VARset]])
  ) |>
  ungroup()

saveRDS(design_grid,
  file = file.path(DATA_DIR, "sim_conditions_singlelevel.rds")
)
message(sprintf("Design grid saved (%d conditions)", nrow(design_grid))) # 192

# -------------------------------------------------------------------------
# 2.  SIMULATION -----------------------------------------------------------
# -------------------------------------------------------------------------
if (RUN_SINGLE_SIM) {
  source("simulate_data.R", local = TRUE) # defines simulate_all_conditions_var1()
  simulate_all_conditions_var1(
    sim_conditions_df = design_grid,
    output_dir        = DATA_DIR
  )
} else {
  message(">> Skipping data generation.")
}

# -------------------------------------------------------------------------
# 3.  QUICK VISUAL CHECKS --------------------------------------------------
# -------------------------------------------------------------------------
if (RUN_SINGLE_CHECK) {
  source("check_simulations.R", local = TRUE)
  run_post_sim_checks_var1(
    data_dir   = DATA_DIR,
    checks_dir = CHECKS_DIR
  )
}

# -------------------------------------------------------------------------
# 4.  STAN FITTING (Normal vs. SkewNormal) --------------------------------
# -------------------------------------------------------------------------
if (RUN_SINGLE_FIT) {
  source("fit_models.R", local = TRUE) # defines fit_var1_copula_models()
  fit_var1_copula_models(
    data_dir            = DATA_DIR,
    fits_dir            = FITS_DIR,
    stan_models_dir     = file.path(BASE_DIR, "stan"),
    sim_conditions_file = file.path(DATA_DIR, "sim_conditions_singlelevel.rds"),
    stan_iter           = 4000,
    stan_warmup         = 2000,
    stan_chains         = 4,
    num_cores           = max(1, NUM_CORES)
  )
}

# -------------------------------------------------------------------------
# 5.  Place‑holders for later phases --------------------------------------
# -------------------------------------------------------------------------
if (RUN_EXP_TEST) {
  # source("run_exponential_test.R")
}
if (RUN_MULTILEVEL) {
  # source("run_multilevel_pipeline.R")
}

message("Pipeline complete.")
