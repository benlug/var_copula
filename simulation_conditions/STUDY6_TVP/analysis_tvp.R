###########################################################################
# analysis_tvp.R  – Study 6: Time-Varying Parameter Copula-VAR
###########################################################################
#
# Extracts and summarizes results from TVP and Constant model fits.
# NOTE: Fits are saved as compact summary objects (not full stanfit) to save space.
#
# Key metrics:
#   - rho trajectory recovery (RMSE, range recovery)
#   - TVP detection (sigma_z posterior)
#   - VAR parameter recovery (bias, coverage)
#
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── folders ------------------------------------------------------------------
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")

dir.create(RES_DIR, FALSE, TRUE)

# ── helpers ------------------------------------------------------------------
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Extract parameter info from compact summary object
extract_from_summary <- function(fit_obj, param_names) {
  if (is.null(fit_obj) || is.null(fit_obj$summary)) return(NULL)
  
  s <- fit_obj$summary
  available <- rownames(s)
  keep <- intersect(param_names, available)
  
  if (length(keep) == 0) return(NULL)
  
  as.data.frame(s[keep, , drop = FALSE], optional = TRUE) |>
    tibble::rownames_to_column("param")
}

# ── parameter sets -----------------------------------------------------------
CORE_PARAMS <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22")
TVP_PARAMS <- c("z0", "sigma_z", "rho_mean", "rho_sd", "rho_range", 
                "rho_start", "rho_end", "rho_change")
CONST_PARAMS <- c("rho")
SIGMA_NG_PARAMS <- c("sigma[1]", "sigma[2]")
SIGMA_EG_PARAMS <- c("sigma_exp[1]", "sigma_exp[2]")

# ── enumerate fits -----------------------------------------------------------
fits <- list.files(FITS_DIR, "^fit_(TVP_NG|TVP_EG|Const_NG|Const_EG)_cond\\d+_rep\\d+\\.rds$",
                   full.names = TRUE)

if (length(fits) == 0) {
  message("No fit files found in ", FITS_DIR)
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
}

message("Found ", length(fits), " fit files")

# ── load simulation metadata -------------------------------------------------
sim_env <- new.env()
load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_env, inherits = FALSE)) {
    p <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid))
    assign(key, safe_read(p), sim_env)
  }
  get(key, sim_env, inherits = FALSE)
}

# ── main analysis loop -------------------------------------------------------
message("Analyzing ", length(fits), " fits...")
rows <- vector("list", length(fits))

for (i in seq_along(fits)) {
  fp <- fits[i]
  m <- str_match(basename(fp), "fit_(TVP_NG|TVP_EG|Const_NG|Const_EG)_cond(\\d+)_rep(\\d+)\\.rds")
  model <- m[2]
  cid <- as.integer(m[3])
  rid <- as.integer(m[4])
  
  is_tvp <- grepl("^TVP", model)
  
  fit_obj <- safe_read(fp)
  
  # Handle compact summary objects (new format) or errors
  if (inherits(fit_obj, "error")) {
    stat <- "bad_rds"
  } else if (is.list(fit_obj) && !is.null(fit_obj$status)) {
    stat <- fit_obj$status
  } else {
    stat <- "unknown_format"
  }
  
  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rows[[i]] <- tibble(model, condition_id = cid, rep_id = rid,
                        param = NA_character_, status = "sim_missing")
    next
  }
  
  # Ground truth
  truth_core <- c(
    `mu[1]` = 0, `mu[2]` = 0,
    phi11 = sim$phi_matrix[1, 1], phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1], phi22 = sim$phi_matrix[2, 2]
  )
  
  # Determine if this is an EG model
  is_eg <- grepl("EG", model)
  
  # Marginal truth (sigma = 1 for standardized data)
  if (is_eg) {
    truth_sigma <- c(`sigma_exp[1]` = 1, `sigma_exp[2]` = 1)
    sigma_params <- SIGMA_EG_PARAMS
  } else {
    truth_sigma <- c(`sigma[1]` = 1, `sigma[2]` = 1)
    sigma_params <- SIGMA_NG_PARAMS
  }
  
  # TVP-specific truth
  if (is_tvp) {
    truth_tvp <- c(
      rho_mean = sim$true_params$rho_mean,
      rho_sd = sim$true_params$rho_sd,
      rho_start = sim$true_params$rho_start,
      rho_end = sim$true_params$rho_end,
      rho_range = sim$true_params$rho_range
    )
  } else {
    # For constant model, rho_mean is the target
    truth_const <- c(rho = sim$true_params$rho_mean)
  }
  
  # Combine truth vectors
  if (is_tvp) {
    truth <- c(truth_core, truth_sigma, truth_tvp)
    want <- c(CORE_PARAMS, sigma_params, TVP_PARAMS)
  } else {
    truth <- c(truth_core, truth_sigma, truth_const)
    want <- c(CORE_PARAMS, sigma_params, CONST_PARAMS)
  }
  
  # Extract summaries from compact object
  if (stat == "ok") {
    sm <- extract_from_summary(fit_obj, want)
  } else {
    sm <- NULL
  }
  
  # Get diagnostics from compact object
  rhat <- if (stat == "ok") fit_obj$max_rhat %||% NA_real_ else NA_real_
  n_div <- if (stat == "ok") fit_obj$n_div %||% NA_integer_ else NA_integer_
  
  rows[[i]] <- if (is.null(sm)) {
    tibble(model, condition_id = cid, rep_id = rid,
           param = NA_character_, status = stat,
           n_div = n_div, max_rhat = rhat)
  } else {
    sm |>
      mutate(
        truth = truth[param],
        bias = ifelse(is.na(truth), NA_real_, mean - truth),
        rel_bias = ifelse(is.na(truth), NA_real_,
                          ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth))),
        cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
        model = model, condition_id = cid, rep_id = rid,
        n_div = n_div, max_rhat = rhat, status = stat
      ) |>
      select(model, condition_id, rep_id, param,
             post_mean = mean, post_sd = sd,
             l95 = `2.5%`, u95 = `97.5%`,
             truth, bias, rel_bias, cover95,
             n_div, max_rhat, status)
  }
  
  if (i %% 200 == 0) message(sprintf("... processed %d/%d", i, length(fits)))
}

rep_tbl <- bind_rows(rows)

if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) {
  rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
}

write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))

message("Fit status summary:")
print(table(Status = rep_tbl$status, Model = rep_tbl$model, useNA = "ifany"))

# ── Condition-level aggregation ----------------------------------------------
good <- rep_tbl$status == "ok" & !is.na(rep_tbl$param)

if (any(good)) {
  # Load design for condition metadata
  design <- readRDS(file.path(DATA_DIR, "sim_conditions.rds"))
  
  cond <- rep_tbl |>
    filter(good) |>
    left_join(
      design |> select(condition_id, T, tvp_pattern, margin_type, direction, VARset),
      by = "condition_id"
    ) |>
    group_by(condition_id, model, param, T, tvp_pattern, margin_type, direction, VARset) |>
    summarise(
      N_valid = n(),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_rel_bias = mean(rel_bias, na.rm = TRUE),
      coverage_95 = mean(cover95, na.rm = TRUE),
      mean_post_sd = mean(post_sd, na.rm = TRUE),
      emp_sd = sd(post_mean, na.rm = TRUE),
      sd_bias = mean_post_sd - emp_sd,
      RMSE = sqrt(mean(bias^2, na.rm = TRUE)),
      mean_n_div = mean(n_div, na.rm = TRUE),
      prop_div = mean(n_div > 0, na.rm = TRUE),
      mean_rhat = mean(max_rhat, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))
  
  write_csv(cond, file.path(RES_DIR, "summary_conditions.csv"))
  message("✓ wrote ", nrow(cond), " condition rows to summary_conditions.csv")
  
  # ── TVP-specific summaries -----------------------------------------------
  # Extract sigma_z estimates for TVP detection analysis
  sigma_z_summary <- rep_tbl |>
    filter(good, param == "sigma_z", grepl("^TVP", model)) |>
    left_join(
      design |> select(condition_id, T, tvp_pattern),
      by = "condition_id"
    ) |>
    group_by(condition_id, model, T, tvp_pattern) |>
    summarise(
      N = n(),
      sigma_z_mean = mean(post_mean, na.rm = TRUE),
      sigma_z_median = median(post_mean, na.rm = TRUE),
      sigma_z_sd = sd(post_mean, na.rm = TRUE),
      prop_sigma_z_gt_0.01 = mean(l95 > 0.01, na.rm = TRUE),
      prop_sigma_z_gt_0.02 = mean(l95 > 0.02, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(sigma_z_summary, file.path(RES_DIR, "summary_sigma_z.csv"))
  message("✓ wrote sigma_z summary to summary_sigma_z.csv")
  
  # ── rho trajectory recovery -----------------------------------------------
  rho_recovery <- rep_tbl |>
    filter(good, param %in% c("rho_mean", "rho_sd", "rho_range", "rho_start", "rho_end")) |>
    left_join(
      design |> select(condition_id, T, tvp_pattern),
      by = "condition_id"
    ) |>
    group_by(condition_id, model, param, T, tvp_pattern) |>
    summarise(
      N = n(),
      mean_estimate = mean(post_mean, na.rm = TRUE),
      mean_truth = mean(truth, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      coverage_95 = mean(cover95, na.rm = TRUE),
      RMSE = sqrt(mean(bias^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  write_csv(rho_recovery, file.path(RES_DIR, "summary_rho_recovery.csv"))
  message("✓ wrote rho recovery summary to summary_rho_recovery.csv")
  
} else {
  message("⚠ no usable draws found in any fit.")
}

message("\n✓ Analysis complete. Results in: ", RES_DIR)
