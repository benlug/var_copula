###########################################################################
# analysis_singlelevel.R  – single‑CPU, robust, analyses alpha
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── folders --------------------------------------------------------------
# Ensure 'this.path' is available or use an alternative way to set BASE_DIR
if (!requireNamespace("this.path", quietly = TRUE)) {
  # Fallback if this.path is not installed
  message("Package 'this.path' not found. Assuming current working directory.")
  BASE_DIR <- getwd()
} else {
  # Use tryCatch in case this.path fails (e.g., when sourced interactively without a file)
  tryCatch(
    {
      BASE_DIR <- this.path::this.dir()
    },
    error = function(e) {
      message("this.path::this.dir() failed. Assuming current working directory.")
      BASE_DIR <<- getwd()
    }
  )
}

DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, FALSE, TRUE)

# ── helpers --------------------------------------------------------------
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Robust helpers with error checking for potentially corrupt/empty fits
count_div <- function(fit) {
  if (!inherits(fit, "stanfit")) {
    return(NA_integer_)
  }
  # Check if sampler params exist before accessing
  sp <- try(get_sampler_params(fit, FALSE), silent = TRUE)
  if (inherits(sp, "try-error") || length(sp) == 0) {
    return(NA_integer_)
  }
  sum(vapply(
    sp,
    \(x) as.integer(sum(x[, "divergent__"])), 0L
  ))
}

max_rhat <- function(fit) {
  if (!inherits(fit, "stanfit")) {
    return(NA_real_)
  }
  # Check if summary exists
  s <- try(summary(fit)$summary, silent = TRUE)
  if (inherits(s, "try-error") || is.null(s) || !("Rhat" %in% colnames(s))) {
    return(NA_real_)
  }
  # Use suppressWarnings as Rhat might be NA if variance is zero
  suppressWarnings(max(s[, "Rhat"], na.rm = TRUE))
}

quick_summary <- function(fit, want_full) {
  base <- unique(sub("\\[.*$", "", want_full))
  # Check against model_pars (includes parameters, transformed parameters, generated quantities)
  keep <- intersect(base, fit@model_pars)

  if (!length(keep)) {
    return(NULL)
  }
  # Use pars argument for efficiency, wrap in try()
  s <- try(summary(fit, pars = keep)$summary, silent = TRUE)
  if (inherits(s, "try-error")) {
    return(NULL)
  }

  s |>
    as.data.frame(optional = TRUE) |>
    tibble::rownames_to_column("param") |>
    filter(param %in% want_full)
}

# ── parameter sets -------------------------------------------------------
CORE <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
SIGMA_NG <- c("sigma[1]", "sigma[2]")
# We analyze the DP parameters (alpha, omega) which should be available
# if calculated in the Stan model's transformed parameters or generated quantities.
EXTRA_SG <- c("omega[1]", "omega[2]", "alpha[1]", "alpha[2]")

# ── enumerate fits -------------------------------------------------------
fits <- list.files(FITS_DIR, "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)

if (length(fits) == 0) {
  message("No fit files found in ", FITS_DIR)
  # Exit if run as a script, otherwise return silently
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

# ── simulation cache -----------------------------------------------------
sim_env <- new.env()

load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_env, inherits = FALSE)) {
    p <- file.path(
      DATA_DIR,
      sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid)
    )
    assign(key, safe_read(p), sim_env)
  }
  get(key, sim_env, inherits = FALSE)
}

message("⏳ analysing ", length(fits), " fits on ONE core …")
rows <- vector("list", length(fits))

for (i in seq_along(fits)) {
  fp <- fits[i]
  m <- str_match(basename(fp), "fit_(SG|NG)_cond(\\d+)_rep(\\d+)\\.rds")
  model <- m[2]
  cid <- as.integer(m[3])
  rid <- as.integer(m[4])

  fit <- safe_read(fp)
  stat <- if (inherits(fit, "stanfit")) {
    # Robust check if sampling actually occurred, iterations > 0, and samples exist
    if (length(fit@sim) > 0 && !is.null(fit@sim$iter) && fit@sim$iter > 0 && length(fit@sim$samples) > 0) {
      "ok"
    } else {
      "empty"
    }
  } else if (inherits(fit, "error")) "bad_rds" else "not_stanfit"

  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rows[[i]] <- tibble(model,
      condition_id = cid, rep_id = rid,
      param = NA_character_, status = "sim_missing"
    )
    next
  }

  # ---------- ground truth ------------------------------------------------
  # FIX: Handle ground truth calculation when the DGP is misspecified (e.g. ChiSq)

  # Helper to make delta from alpha, handling NA
  make_delta <- function(a) {
    if (is.na(a)) {
      return(NA_real_)
    }
    a / sqrt(1 + a^2)
  }

  # Helper to safely extract alpha ONLY if the type is skewnormal
  get_sn_alpha <- function(margin_params) {
    # Check if the simulation used skewnormal distribution
    if (is.null(margin_params) || (!is.null(margin_params$type) && margin_params$type != "skewnormal")) {
      return(NA_real_) # Not applicable if not skew-normal (e.g. ChiSq)
    }
    # Return alpha if present, otherwise NA (shouldn't happen if type is skewnormal)
    return(margin_params$alpha %||% NA_real_)
  }

  alpha1 <- get_sn_alpha(sim$true_params$margin1)
  alpha2 <- get_sn_alpha(sim$true_params$margin2)

  # Function to calculate omega (scale for standardized variance=1), handling NA delta
  calc_omega <- function(d) {
    if (is.na(d)) {
      return(NA_real_)
    }
    # Formula for omega when variance is 1: sqrt(1 / (1 - 2*delta^2/pi))
    denom <- 1 - 2 * d^2 / pi
    # Safety check for denominator near zero
    if (denom <= 1e-9) {
      return(NA_real_)
    }
    sqrt(1 / denom)
  }

  delta1 <- make_delta(alpha1)
  delta2 <- make_delta(alpha2)
  omega1 <- calc_omega(delta1)
  omega2 <- calc_omega(delta2)

  # Define the truth vector
  truth <- c(
    # Since simulations are standardized (Mean=0), true mu is 0.
    `mu[1]` = 0, `mu[2]` = 0,
    phi11 = sim$phi_matrix[1, 1], phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1], phi22 = sim$phi_matrix[2, 2],
    rho = sim$rho
  )

  if (model == "NG") {
    # Since simulations are standardized (Var=1), true sigma is 1 for the NG model.
    truth <- c(truth, `sigma[1]` = 1, `sigma[2]` = 1)
  }

  if (model == "SG") {
    # For SG, we use the calculated (or NA) values for the DP parameters.
    truth <- c(truth,
      `omega[1]` = omega1, `omega[2]` = omega2,
      `alpha[1]` = alpha1, `alpha[2]` = alpha2
    )
  }

  want <- c(
    CORE,
    if (model == "NG") SIGMA_NG,
    if (model == "SG") EXTRA_SG
  )

  sm <- if (stat == "ok") quick_summary(fit, want) else NULL
  rhat <- if (stat == "ok") max_rhat(fit) else NA_real_
  n_div <- count_div(fit) # Calculate divergences if possible

  rows[[i]] <- if (is.null(sm)) {
    tibble(model,
      condition_id = cid, rep_id = rid,
      param = NA_character_, status = stat,
      n_div = n_div, max_rhat = rhat
    )
  } else {
    sm |>
      mutate(
        truth = truth[param],
        # FIX: Handle NA in truth for bias calculation
        bias = ifelse(is.na(truth), NA_real_, mean - truth),

        # FIX: Handle relative bias calculation carefully when truth is 0 (e.g., mu) or NA
        rel_bias = ifelse(is.na(truth), NA_real_,
          ifelse(abs(truth) < .Machine$double.eps,
            bias, bias / abs(truth)
          )
        ),

        # FIX: Handle NA in truth for coverage calculation
        cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
        model = model, condition_id = cid, rep_id = rid,
        n_div = n_div,
        max_rhat = rhat,
        status = stat
      ) |>
      # Ensure columns are in a consistent order
      select(model, condition_id, rep_id, param,
        post_mean = mean, post_sd = sd,
        l95 = `2.5%`, u95 = `97.5%`,
        truth, bias, rel_bias, cover95,
        n_div, max_rhat, status
      )
  }

  if (i %% 500 == 0) message(sprintf("... processed %d/%d", i, length(fits)))
}

rep_tbl <- bind_rows(rows)

# Ensure logical/numeric types are correct after binding
# (especially for cover95 which might become logical if all NA or only TRUE/FALSE)
if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) {
  # Convert logical (TRUE/FALSE/NA) to numeric (1/0/NA) for aggregation (mean)
  rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
}


write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))

message("Fit status summary:")
print(table(rep_tbl$status, useNA = "ifany"))

# Aggregate results only for successfully fitted models with parameters
good <- rep_tbl$status == "ok" & !is.na(rep_tbl$param)
if (any(good)) {
  cond <- rep_tbl |>
    filter(good) |>
    group_by(condition_id, model, param) |>
    summarise(
      N_valid = n(), # Count of valid replications
      N_truth_avail = sum(!is.na(truth)), # Count where truth is available
      # Use na.rm=TRUE for aggregations where truth might be NA
      mean_bias = mean(bias, na.rm = TRUE),
      mean_rel_bias = mean(rel_bias, na.rm = TRUE),
      coverage_95 = mean(cover95, na.rm = TRUE),
      mean_post_sd = mean(post_sd, na.rm = TRUE),
      # Empirical SD calculation uses post_mean
      emp_sd = sd(post_mean, na.rm = TRUE),
      # SD-Bias: difference between average estimated SD and empirical SD
      sd_bias = mean_post_sd - emp_sd,
      mean_n_div = mean(n_div, na.rm = TRUE),
      prop_div = mean(n_div > 0, na.rm = TRUE), # Proportion of runs with divergences
      mean_rhat = mean(max_rhat, na.rm = TRUE),
      .groups = "drop"
    )

  # Replace NaN resulting from mean(..., na.rm=TRUE) on all-NA inputs (e.g. if N_truth_avail=0) with NA_real_
  cond <- cond |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

  write_csv(cond, file.path(RES_DIR, "summary_conditions.csv"))
  message("✓ wrote ", nrow(cond), " condition rows to summary_conditions.csv")
} else {
  message("⚠ no usable draws found in any fit.")
}
