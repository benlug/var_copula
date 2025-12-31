###########################################################################
# analysis_singlelevel.R  – Analysis of compact fit summaries
#   Matches TVP/Study 2 format: fit_summary$summary is a matrix, not data.frame
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── Folders ---------------------------------------------------------------
if (!requireNamespace("this.path", quietly = TRUE)) {
  BASE_DIR <- getwd()
} else {
  tryCatch({
    BASE_DIR <- this.path::this.dir()
  }, error = function(e) {
    BASE_DIR <<- getwd()
  })
}

DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR  <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, FALSE, TRUE)

# ── Helpers ---------------------------------------------------------------
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Parameter sets --------------------------------------------------------
CORE <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
SIGMA_NG <- c("sigma[1]", "sigma[2]")
EXTRA_SG <- c("omega[1]", "omega[2]", "alpha[1]", "alpha[2]")

# ── Enumerate expected fits (include missing) -----------------------------
# Goal: avoid selection bias from silently excluding failed/missing fits.
design_path <- file.path(DATA_DIR, "sim_conditions_singlelevel.rds")
design <- if (file.exists(design_path)) safe_read(design_path) else NULL

expected <- NULL

if (!is.null(design) && !inherits(design, "error") &&
    all(c("condition_id", "n_reps") %in% names(design))) {

  expected_base <- design |>
    transmute(
      condition_id = as.integer(condition_id),
      n_reps = as.integer(n_reps)
    ) |>
    tidyr::uncount(weights = n_reps, .id = "rep_id") |>
    select(condition_id, rep_id)

  expected <- tidyr::crossing(expected_base, model = c("SG", "NG")) |>
    arrange(condition_id, rep_id, model) |>
    mutate(
      fit_path = file.path(
        FITS_DIR,
        sprintf("fit_%s_cond%03d_rep%03d.rds", model, condition_id, rep_id)
      )
    )
} else {
  # Fallback: only enumerate existing fit files
  fits <- list.files(FITS_DIR, "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (length(fits) == 0) {
    message("No fit files found in ", FITS_DIR)
    if (sys.nframe() == 0) q(status = 1) else return(invisible())
  }
  meta <- str_match(basename(fits), "fit_(SG|NG)_cond(\\d+)_rep(\\d+)\\.rds")
  expected <- tibble(
    model = meta[, 2],
    condition_id = as.integer(meta[, 3]),
    rep_id = as.integer(meta[, 4]),
    fit_path = fits
  ) |>
    arrange(condition_id, rep_id, model)
}

# ── Simulation cache ------------------------------------------------------
sim_env <- new.env()

load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_env, inherits = FALSE)) {
    p <- file.path(DATA_DIR,
                   sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid))
    assign(key, safe_read(p), sim_env)
  }
  get(key, sim_env, inherits = FALSE)
}

message("Analysing ", nrow(expected), " expected fits (including missing)...")
rows <- vector("list", nrow(expected))

for (i in seq_len(nrow(expected))) {
  fp <- expected$fit_path[i]
  model <- as.character(expected$model[i])
  cid <- as.integer(expected$condition_id[i])
  rid <- as.integer(expected$rep_id[i])

  # Missing fit file
  if (!file.exists(fp)) {
    # Backwards-compatibility: older fitting scripts wrote failures to
    # fit_ERR_* / fit_MISC_* without creating the canonical fit_* file.
    fp_err  <- file.path(FITS_DIR, sprintf("fit_ERR_%s_cond%03d_rep%03d.rds", model, cid, rid))
    fp_misc <- file.path(FITS_DIR, sprintf("fit_MISC_%s_cond%03d_rep%03d.rds", model, cid, rid))
    legacy_status <- if (file.exists(fp_err)) "legacy_error" else if (file.exists(fp_misc)) "legacy_misc" else "missing_fit"

    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = legacy_status,
      n_div = NA_integer_, max_rhat = NA_real_
    )
    next
  }

  fit_data <- safe_read(fp)
  
  # Handle read errors
  if (inherits(fit_data, "error")) {
    rows[[i]] <- tibble(model, condition_id = cid, rep_id = rid,
                        param = NA_character_, status = "bad_rds")
    next
  }
  
  # Extract status
  stat <- fit_data$status %||% "unknown"
  
  # Load simulation data for ground truth
  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rows[[i]] <- tibble(model, condition_id = cid, rep_id = rid,
                        param = NA_character_, status = "sim_missing")
    next
  }

  # ── Ground truth calculation --------------------------------------------
  make_delta <- function(a) {
    if (is.na(a)) return(NA_real_)
    a / sqrt(1 + a^2)
  }

  get_sn_alpha <- function(margin_params) {
    if (is.null(margin_params) || 
        (!is.null(margin_params$type) && margin_params$type != "skewnormal")) {
      return(NA_real_)
    }
    return(margin_params$alpha %||% NA_real_)
  }

  alpha1 <- get_sn_alpha(sim$true_params$margin1)
  alpha2 <- get_sn_alpha(sim$true_params$margin2)

  calc_omega <- function(d) {
    if (is.na(d)) return(NA_real_)
    denom <- 1 - 2 * d^2 / pi
    if (denom <= 1e-9) return(NA_real_)
    sqrt(1 / denom)
  }

  delta1 <- make_delta(alpha1)
  delta2 <- make_delta(alpha2)
  omega1 <- calc_omega(delta1)
  omega2 <- calc_omega(delta2)

  # Build truth vector
  truth <- c(
    `mu[1]` = 0, `mu[2]` = 0,
    phi11 = sim$phi_matrix[1, 1], phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1], phi22 = sim$phi_matrix[2, 2],
    rho = sim$rho
  )

  if (model == "NG") {
    truth <- c(truth, `sigma[1]` = 1, `sigma[2]` = 1)
  }

  if (model == "SG") {
    truth <- c(truth,
               `omega[1]` = omega1, `omega[2]` = omega2,
               `alpha[1]` = alpha1, `alpha[2]` = alpha2)
  }

  want <- c(CORE,
            if (model == "NG") SIGMA_NG,
            if (model == "SG") EXTRA_SG)

  # Extract summary - it's a matrix in TVP format, not a data frame
  sm_raw <- fit_data$summary
  n_div <- fit_data$n_div %||% NA_integer_
  max_rhat <- fit_data$max_rhat %||% NA_real_
  
  if (is.null(sm_raw) || stat != "ok") {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = stat,
      n_div = n_div, max_rhat = max_rhat
    )
    next
  }

  # Convert matrix to data frame
  sm <- as.data.frame(sm_raw)
  sm$param <- rownames(sm_raw)
  rownames(sm) <- NULL
  
  # Filter to wanted parameters
  sm <- sm[sm$param %in% want, , drop = FALSE]
  
  if (nrow(sm) == 0) {
    rows[[i]] <- tibble(model, condition_id = cid, rep_id = rid,
                        param = NA_character_, status = "no_params")
    next
  }

  # Build result
  result <- data.frame(
    model = model,
    condition_id = cid,
    rep_id = rid,
    param = sm$param,
    post_mean = sm$mean,
    post_sd = sm$sd,
    l95 = sm$`2.5%`,
    u95 = sm$`97.5%`,
    stringsAsFactors = FALSE
  )
  
  result$truth <- truth[result$param]
  result$bias <- ifelse(is.na(result$truth), NA_real_, result$post_mean - result$truth)
  result$rel_bias <- ifelse(is.na(result$truth), NA_real_,
                            ifelse(abs(result$truth) < .Machine$double.eps,
                                   result$bias, result$bias / abs(result$truth)))
  result$cover95 <- ifelse(is.na(result$truth), NA,
                           result$l95 <= result$truth & result$u95 >= result$truth)
  result$n_div <- n_div
  result$max_rhat <- max_rhat
  result$status <- stat
  
  rows[[i]] <- result

  if (i %% 500 == 0) message(sprintf("... processed %d/%d", i, nrow(expected)))
}

rep_tbl <- bind_rows(rows)

# Ensure proper types
if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) {
  rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
}

write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))

# One-row-per-fit status table (includes missing fits)
fit_status <- rep_tbl |>
  distinct(model, condition_id, rep_id, status, n_div, max_rhat)

write_csv(fit_status, file.path(RES_DIR, "fit_status_by_replication.csv"))

fit_status_cond <- fit_status |>
  group_by(condition_id, model) |>
  summarise(
    N_expected = n(),
    N_ok = sum(status == "ok"),
    N_missing_fit = sum(status == "missing_fit"),
    N_not_ok = sum(status != "ok" & status != "missing_fit"),
    prop_ok = ifelse(N_expected > 0, N_ok / N_expected, NA_real_),
    .groups = "drop"
  )

write_csv(fit_status_cond, file.path(RES_DIR, "fit_status_by_condition.csv"))

message("Fit status summary:")
print(table(rep_tbl$status, useNA = "ifany"))

# Aggregate results
good <- rep_tbl$status == "ok" & !is.na(rep_tbl$param)
if (any(good)) {
  cond <- rep_tbl |>
    filter(good) |>
    group_by(condition_id, model, param) |>
    summarise(
      N_valid = n(),
      N_truth_avail = sum(!is.na(truth)),
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

  # Replace NaN with NA
  cond <- cond |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

  # Attach denominators for transparency (same for all params within a cell)
  cond <- cond |>
    left_join(fit_status_cond, by = c("condition_id", "model"))

  write_csv(cond, file.path(RES_DIR, "summary_conditions.csv"))
  message("Wrote ", nrow(cond), " condition rows to summary_conditions.csv")
} else {
  message("No usable draws found in any fit.")
}
