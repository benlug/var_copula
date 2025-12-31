###########################################################################
# analysis_singlelevel.R  – (Updated for Study 2: Includes EG)
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── folders --------------------------------------------------------------
# Robust directory detection (assuming this script is run from the project root)
BASE_DIR <- getwd()

# Check if Study 2 directories exist (created by the updated run_pipeline.R), otherwise fall back to Study 1
if (dir.exists(file.path(BASE_DIR, "data"))) {
  message("Analyzing Study 2 results.")
  DATA_DIR <- file.path(BASE_DIR, "data")
  FITS_DIR <- file.path(BASE_DIR, "fits")
  RES_DIR <- file.path(BASE_DIR, "results")
} else {
  message("Analyzing Study 1 results.")
  DATA_DIR <- file.path(BASE_DIR, "data")
  FITS_DIR <- file.path(BASE_DIR, "fits")
  RES_DIR <- file.path(BASE_DIR, "results")
}

dir.create(RES_DIR, FALSE, TRUE)

# ── helpers --------------------------------------------------------------
# (Helpers remain the same as the corrected Study 1 version)
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

count_div <- function(fit) {
  # Lightweight fit artifact (saved by fit_models.R in "summary" mode)
  if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    return(as.integer(fit$n_div %||% NA_integer_))
  }
  if (!inherits(fit, "stanfit")) {
    return(NA_integer_)
  }
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
  # Lightweight fit artifact
  if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    s <- fit$summary
    if (is.data.frame(s) && ("Rhat" %in% names(s))) {
      r <- suppressWarnings(as.numeric(s$Rhat))
      r <- r[is.finite(r)]
      return(if (length(r)) max(r) else NA_real_)
    }
    return(NA_real_)
  }
  if (!inherits(fit, "stanfit")) {
    return(NA_real_)
  }
  s <- try(summary(fit)$summary, silent = TRUE)
  if (inherits(s, "try-error") || is.null(s) || !("Rhat" %in% colnames(s))) {
    return(NA_real_)
  }
  suppressWarnings(max(s[, "Rhat"], na.rm = TRUE))
}

quick_summary <- function(fit, want_full) {
  # Lightweight fit artifact
  if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    s <- fit$summary
    if (!is.data.frame(s)) return(NULL)
    if (!"param" %in% names(s)) {
      s <- tibble::rownames_to_column(s, "param")
    }
    return(dplyr::filter(s, param %in% want_full))
  }
  base <- unique(sub("\\[.*$", "", want_full))
  keep <- intersect(base, fit@model_pars)
  if (!length(keep)) {
    return(NULL)
  }
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
EXTRA_SG <- c("omega[1]", "omega[2]", "alpha[1]", "alpha[2]")
# STUDY 2 ADDITION: Exponential parameters (Scale)
EXTRA_EG <- c("sigma_exp[1]", "sigma_exp[2]")

# ── enumerate fits -------------------------------------------------------
# Update regex to include EG
fits <- list.files(FITS_DIR, "^fit_(SG|NG|EG)_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)

if (length(fits) == 0) {
  message("No fit files found in ", FITS_DIR)
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

# ── simulation cache -----------------------------------------------------
# (load_sim remains the same)
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
  # Update regex to include EG
  m <- str_match(basename(fp), "fit_(SG|NG|EG)_cond(\\d+)_rep(\\d+)\\.rds")
  model <- m[2]
  cid <- as.integer(m[3])
  rid <- as.integer(m[4])

  fit <- safe_read(fp)
  # Fit artifacts can be either full stanfit objects or lightweight summaries
  stat <- if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    if (is.data.frame(fit$summary) && nrow(fit$summary) > 0) "ok" else "empty"
  } else if (inherits(fit, "stanfit")) {
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
  # We use the logic established in the Study 1 review to handle misspecification.

  # Helpers for Skew-Normal (will return NA if DGP is not Skew-Normal)
  make_delta <- function(a) {
    if (is.na(a)) {
      return(NA_real_)
    }
    a / sqrt(1 + a^2)
  }

  get_sn_alpha <- function(margin_params) {
    if (is.null(margin_params) || (!is.null(margin_params$type) && margin_params$type != "skewnormal")) {
      return(NA_real_)
    }
    return(margin_params$alpha %||% NA_real_)
  }

  calc_omega <- function(d) {
    if (is.na(d)) {
      return(NA_real_)
    }
    denom <- 1 - 2 * d^2 / pi
    if (denom <= 1e-9) {
      return(NA_real_)
    }
    sqrt(1 / denom)
  }

  alpha1 <- get_sn_alpha(sim$true_params$margin1)
  alpha2 <- get_sn_alpha(sim$true_params$margin2)
  omega1 <- calc_omega(make_delta(alpha1))
  omega2 <- calc_omega(make_delta(alpha2))

  # STUDY 2: Exponential Ground Truth
  # Check if the DGP was exponential. If so, the true standardized scale (sigma_exp) is 1.
  is_exponential_dgp <- (!is.null(sim$true_params$margin1$type) && sim$true_params$margin1$type == "exponential")
  true_sigma_exp <- if (is_exponential_dgp) 1 else NA_real_

  # Define the truth vector (CORE parameters)
  # --- Effective copula rho (accounts for mirroring) --------------------
  # Mirroring exactly one margin is a monotone-decreasing transform of its copula U,
  # which flips the sign of the Gaussian copula parameter on the latent-normal scale.
  # Target for evaluation is rho_eff = s1 * s2 * rho, where s_j = -1 if mirrored.
  m1_mirror <- sim$true_params$margin1$mirror %||% FALSE
  m2_mirror <- sim$true_params$margin2$mirror %||% FALSE
  s1 <- if (isTRUE(m1_mirror)) -1 else 1
  s2 <- if (isTRUE(m2_mirror)) -1 else 1
  rho_eff <- s1 * s2 * sim$rho

  # Define the truth vector (CORE parameters)
  truth <- c(
    `mu[1]` = 0, `mu[2]` = 0, # Standardized data
    phi11 = sim$phi_matrix[1, 1], phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1], phi22 = sim$phi_matrix[2, 2],
    rho = rho_eff
  )

  # Append model-specific parameters
  if (model == "NG") {
    # Standardized data (Var=1), so true sigma is 1.
    truth <- c(truth, `sigma[1]` = 1, `sigma[2]` = 1)
  }

  if (model == "SG") {
    # If DGP is Exponential, these will be NA (correctly handling misspecification)
    truth <- c(truth,
      `omega[1]` = omega1, `omega[2]` = omega2,
      `alpha[1]` = alpha1, `alpha[2]` = alpha2
    )
  }

  if (model == "EG") {
    # STUDY 2 ADDITION
    truth <- c(truth, `sigma_exp[1]` = true_sigma_exp, `sigma_exp[2]` = true_sigma_exp)
  }

  want <- c(
    CORE,
    if (model == "NG") SIGMA_NG,
    if (model == "SG") EXTRA_SG,
    if (model == "EG") EXTRA_EG
  ) # STUDY 2 ADDITION

  sm <- if (stat == "ok") quick_summary(fit, want) else NULL
  rhat <- if (stat == "ok") max_rhat(fit) else NA_real_
  n_div <- count_div(fit)

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
        # Bias and coverage calculations (already robust to NA truth)
        bias = ifelse(is.na(truth), NA_real_, mean - truth),
        rel_bias = ifelse(is.na(truth), NA_real_,
          ifelse(abs(truth) < .Machine$double.eps,
            bias, bias / abs(truth)
          )
        ),
        cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
        model = model, condition_id = cid, rep_id = rid,
        n_div = n_div,
        max_rhat = rhat,
        status = stat
      ) |>
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

# Ensure cover95 is numeric for aggregation
if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) {
  rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
}


write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))

message("Fit status summary:")
print(table(Status = rep_tbl$status, Model = rep_tbl$model, useNA = "ifany"))

# Aggregate results
good <- rep_tbl$status == "ok" & !is.na(rep_tbl$param)
if (any(good)) {
  first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) NA_real_ else x[1]
  }

  cond <- rep_tbl |>
    filter(good) |>
    group_by(condition_id, model, param) |>
    summarise(
      N_valid = n(),
      N_truth_avail = sum(!is.na(truth)),
      truth = first_non_na(truth),
      # Use na.rm=TRUE as truth might be NA for misspecified models
      mean_bias = mean(bias, na.rm = TRUE),
      mean_rel_bias = mean(rel_bias, na.rm = TRUE),
      coverage_95 = mean(cover95, na.rm = TRUE),
      mean_post_sd = mean(post_sd, na.rm = TRUE),
      emp_sd = sd(post_mean, na.rm = TRUE),
      sd_bias = mean_post_sd - emp_sd,
      mean_n_div = mean(n_div, na.rm = TRUE),
      prop_div = mean(n_div > 0, na.rm = TRUE),
      mean_rhat = mean(max_rhat, na.rm = TRUE),
      .groups = "drop"
    )

  # Clean up potential NaN results
  cond <- cond |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

  write_csv(cond, file.path(RES_DIR, "summary_conditions.csv"))
  message("✓ wrote ", nrow(cond), " condition rows to summary_conditions.csv")
} else {
  message("⚠ no usable draws found in any fit.")
}
