###########################################################################
# process_fits.R  – updated 2025‑06‑12
#
# Reads Stan fit objects produced by fit_models.R, merges them with
# simulation metadata, and computes bias / coverage / diagnostics.
#  ▸ NEW: true σ is now available even when the margin is skew‑normal.
###########################################################################

library(rstan)
library(tidyverse)
library(stringr)
library(posterior)
library(this.path)

cat("--- Starting Fit Processing ---\n")

## ------------------------------------------------------------------ ##
##  0.  Paths & meta‑data                                             ##
## ------------------------------------------------------------------ ##
RESULTS_DIR <- tryCatch(this.dir(), error = function(e) getwd())
setwd(RESULTS_DIR)

FITS_DIR <- file.path(RESULTS_DIR, "../fits")
DATA_DIR <- file.path(RESULTS_DIR, "../data")
if (!dir.exists(FITS_DIR)) stop("Fits directory not found: ", FITS_DIR)
if (!dir.exists(DATA_DIR)) stop("Data directory not found: ", DATA_DIR)

SIM_COND_FILE <- file.path(DATA_DIR, "sim_conditions.rds")
if (!file.exists(SIM_COND_FILE)) stop("Simulation conditions file missing.")
sim_conditions_df <- readRDS(SIM_COND_FILE) %>%
  mutate(condition_id = as.integer(condition_id))
cat("Loaded", nrow(sim_conditions_df), "simulation conditions.\n")

`%||%` <- function(a, b) if (!is.null(a)) a else b # null‑coalescing


## ------------------------------------------------------------------ ##
##  1.  Helper functions                                              ##
## ------------------------------------------------------------------ ##

# Safe RDS reader that distinguishes “file missing” / “not stanfit” etc.
safe_read_stanfit <- function(path) {
  if (!file.exists(path)) {
    return(list(fit = NULL, reason = "file missing"))
  }
  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) {
    return(list(fit = NULL, reason = obj$message))
  }
  if (!inherits(obj, "stanfit")) {
    return(list(fit = NULL, reason = "object not stanfit"))
  }
  if (length(obj@sim$samples) == 0) {
    return(list(fit = NULL, reason = "empty stanfit"))
  }
  list(fit = obj, reason = NULL)
}

# Energy‑fraction‑of‑missing‑information (E‑FMI)
calc_e_fmi <- function(energy) {
  energy <- energy[is.finite(energy)]
  n <- length(energy)
  if (n < 3) {
    return(NA_real_)
  }
  mu <- mean(energy)
  v <- var(energy)
  rho1 <- sum((energy[-1] - mu) * (energy[-n] - mu)) / (n - 1)
  denom <- v + 2 * rho1
  ifelse(denom > 0, v / denom, NA_real_)
}

# ------------  TRUE‑VALUE look‑up (patched)  --------------------
get_true_value_var1 <- function(par, tp) {
  # tp = true_params list stored with each simulated data set
  m1dist <- tp$margin1_dist
  m2dist <- tp$margin2_dist
  m1 <- tp$margin1_params
  m2 <- tp$margin2_params

  pop_sd_skewnormal <- function(omega, alpha) {
    delta <- alpha / sqrt(1 + alpha^2)
    omega * sqrt(1 - 2 * delta^2 / pi)
  }

  switch(par,
    "phi11" = tp$phi[1, 1],
    "phi12" = tp$phi[1, 2],
    "phi21" = tp$phi[2, 1],
    "phi22" = tp$phi[2, 2],
    "mu[1]" = tp$mu[1],
    "mu[2]" = tp$mu[2],

    ## ---- residual & shape parameters ----
    "sigma[1]" = {
      if (m1dist == "normal") {
        m1$sd
      } else if (m1dist == "skewnormal") {
        pop_sd_skewnormal(m1$omega, m1$alpha)
      } else {
        NA_real_
      }
    },
    "sigma[2]" = {
      if (m2dist == "normal") {
        m2$sd
      } else if (m2dist == "skewnormal") {
        pop_sd_skewnormal(m2$omega, m2$alpha)
      } else {
        NA_real_
      }
    },
    "xi[1]" = if (m1dist == "skewnormal") m1$xi else NA_real_,
    "xi[2]" = if (m2dist == "skewnormal") m2$xi else NA_real_,
    "omega[1]" = if (m1dist == "skewnormal") m1$omega else NA_real_,
    "omega[2]" = if (m2dist == "skewnormal") m2$omega else NA_real_,
    "alpha[1]" = if (m1dist == "skewnormal") m1$alpha else NA_real_,
    "alpha[2]" = if (m2dist == "skewnormal") m2$alpha else NA_real_,

    ## ---- copula parameters ----
    "rho" = if (tp$copula_type == "gaussian") tp$copula_param_value else NA_real_,
    "theta" = if (tp$copula_type == "clayton") tp$copula_param_value else NA_real_,
    "tau" = tp$copula_tau,
    NA_real_
  )
}

get_param_category_var1 <- function(par) {
  if (grepl("^phi\\d{2}$", par)) {
    "VAR Coeffs"
  } else if (grepl("^mu\\[\\d\\]$", par)) {
    "Intercepts"
  } else if (grepl("^(sigma|xi|omega|alpha)\\[", par)) {
    "Residual Params"
  } else if (par %in% c("rho", "theta", "tau")) {
    "Copula Params"
  } else {
    "Other"
  }
}

## ------------------------------------------------------------------ ##
##  2.  Iterate over all fit files                                    ##
## ------------------------------------------------------------------ ##

fit_files <- list.files(FITS_DIR,
  pattern = "^fit_[A-Z]{2}_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)
if (length(fit_files) == 0) stop("No fit files found.")

param_results <- list()
sampler_results <- list()
fail_log <- list()

for (ff in fit_files) {
  base <- basename(ff)
  mt <- str_match(base, "^fit_([A-Z]{2})_cond(\\d+)_rep(\\d+)\\.rds$")
  mdl <- mt[2]
  cond <- as.integer(mt[3])
  rep <- as.integer(mt[4])

  rd <- safe_read_stanfit(ff)
  if (is.null(rd$fit)) {
    fail_log[[length(fail_log) + 1]] <- c(file = base, reason = rd$reason)
    next
  }
  fit <- rd$fit

  sim_file <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cond, rep))
  if (!file.exists(sim_file)) {
    fail_log[[length(fail_log) + 1]] <- c(file = base, reason = "simulation data missing")
    next
  }
  sim_dat <- readRDS(sim_file)
  true_param <- sim_dat$true_params

  ## ---- posterior summaries ----
  summ <- summary(fit)$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    as_tibble() %>%
    select(parameter,
      post_mean = mean, post_sd = sd,
      ci_low = `2.5%`, ci_high = `97.5%`, n_eff, Rhat
    )

  summ <- summ %>%
    mutate(
      condition_id = cond,
      rep_i = rep,
      fitted_model_code = mdl,
      true_value = map_dbl(parameter, get_true_value_var1, tp = true_param),
      bias = post_mean - true_value,
      coverage = ifelse(!is.na(true_value),
        ci_low <= true_value & true_value <= ci_high, NA
      ),
      rel_bias = ifelse(!is.na(true_value) & abs(true_value) > 1e-6,
        bias / abs(true_value), NA
      ),
      param_category = map_chr(parameter, get_param_category_var1)
    )
  param_results[[length(param_results) + 1]] <- summ

  ## ---- sampler diagnostics ----
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  div_tot <- sum(map_dbl(sp, ~ sum(.x[, "divergent__"])))
  max_td <- attr(sp[[1]], "max_treedepth") %||% 10
  td_hits <- sum(map_dbl(sp, ~ sum(.x[, "treedepth__"] >= max_td)))
  efmi <- mean(map_dbl(sp, ~ calc_e_fmi(.x[, "energy__"])), na.rm = TRUE)
  acc_stat <- mean(map_dbl(sp, ~ mean(.x[, "accept_stat__"])), na.rm = TRUE)

  sampler_results[[length(sampler_results) + 1]] <- tibble(
    condition_id = cond, rep_i = rep, fitted_model_code = mdl,
    divergences = div_tot, maxdepth_exceeded = td_hits,
    eFMI = efmi, avg_accept_stat = acc_stat
  )
}

cat(
  "Finished reading fits: ",
  length(param_results), " successful, ",
  length(fail_log), " failed.\n"
)

## ------------------------------------------------------------------ ##
##  3.  Save aggregated outputs                                       ##
## ------------------------------------------------------------------ ##
if (length(param_results) > 0) {
  param_df <- bind_rows(param_results) %>%
    left_join(
      select(
        sim_conditions_df, condition_id, dgp_copula_type,
        dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22
      ),
      by = "condition_id"
    )
  saveRDS(param_df, file.path(RESULTS_DIR, "parameter_summary.rds"))
  cat("  ▸ parameter_summary.rds written\n")
}

if (length(sampler_results) > 0) {
  samp_df <- bind_rows(sampler_results) %>%
    left_join(
      select(
        sim_conditions_df, condition_id, dgp_copula_type,
        dgp_alpha1, dgp_alpha2, dgp_tau, T, phi11:phi22
      ),
      by = "condition_id"
    )
  saveRDS(samp_df, file.path(RESULTS_DIR, "sampler_summary.rds"))
  cat("  ▸ sampler_summary.rds written\n")
}

if (length(fail_log) > 0) {
  write.csv(bind_rows(fail_log),
    file.path(RESULTS_DIR, "failed_fit_report.csv"),
    row.names = FALSE
  )
  cat(
    "  ▸ failed_fit_report.csv written (",
    length(fail_log), " entries)\n"
  )
}

cat("--- Fit Processing finished ---\n")
