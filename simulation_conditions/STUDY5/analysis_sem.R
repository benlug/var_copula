###########################################################################
# analysis_sem.R — Study 5 (SEM): Analysis of compact fit summaries
#
# Aligns with analysis_singlelevel.R (Studies 1–3):
#   - Reads compact fit summaries produced by fit_models_SEM.R
#   - Enumerates expected fits from the design grid so missing/failed fits
#     are represented (avoids silent selection bias)
#   - Writes:
#       results_sem/summary_replications_sem.csv
#       results_sem/fit_status_by_replication_sem.csv
#       results_sem/fit_status_by_condition_sem.csv
#       results_sem/summary_conditions_sem.csv
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ---- folders ------------------------------------------------------------
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
FITS_DIR <- file.path(BASE_DIR, "fits_sem")
RES_DIR  <- file.path(BASE_DIR, "results_sem")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- parameter sets -----------------------------------------------------
CORE  <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
EXTRA <- c("sigma_exp[1]", "sigma_exp[2]")

# ---- enumerate expected fits (include missing) --------------------------
design_path <- file.path(DATA_DIR, "sim_conditions_sem.rds")
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

  expected <- tidyr::crossing(expected_base, model = c("EI", "EL")) |>
    arrange(condition_id, rep_id, model) |>
    mutate(
      fit_path = file.path(
        FITS_DIR,
        sprintf("fit_%s_cond%03d_rep%03d.rds", model, condition_id, rep_id)
      )
    )

} else {
  # Fallback: only enumerate existing fit files
  fits <- list.files(FITS_DIR, "^fit_(EI|EL)_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(fits)) stop("No SEM fit files found in ", FITS_DIR)
  meta <- stringr::str_match(basename(fits), "fit_(EI|EL)_cond(\\d+)_rep(\\d+)\\.rds")
  expected <- tibble(
    model = meta[, 2],
    condition_id = as.integer(meta[, 3]),
    rep_id = as.integer(meta[, 4]),
    fit_path = fits
  ) |>
    arrange(condition_id, rep_id, model)
}

# ---- simulation cache ---------------------------------------------------
sim_env <- new.env(parent = emptyenv())
load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_env, inherits = FALSE)) {
    p <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid))
    assign(key, safe_read(p), sim_env)
  }
  get(key, sim_env, inherits = FALSE)
}

message("Analysing ", nrow(expected), " expected SEM fits (including missing)…")
rows <- vector("list", nrow(expected))

for (i in seq_len(nrow(expected))) {
  fp <- expected$fit_path[i]
  model <- as.character(expected$model[i])
  cid <- as.integer(expected$condition_id[i])
  rid <- as.integer(expected$rep_id[i])

  # Missing fit file
  if (!file.exists(fp)) {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = "missing_fit",
      n_div = NA_integer_, max_rhat = NA_real_
    )
    next
  }

  fit_data <- safe_read(fp)

  # Handle read errors
  if (inherits(fit_data, "error")) {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = "bad_rds",
      n_div = NA_integer_, max_rhat = NA_real_
    )
    next
  }

  stat <- fit_data$status %||% "unknown"
  n_div <- fit_data$n_div %||% NA_integer_
  max_rhat <- fit_data$max_rhat %||% NA_real_

  # Load simulation for ground truth
  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = "sim_missing",
      n_div = n_div, max_rhat = max_rhat
    )
    next
  }

  truth <- c(
    `mu[1]` = 0, `mu[2]` = 0,
    phi11 = sim$phi_matrix[1, 1],
    phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1],
    phi22 = sim$phi_matrix[2, 2],
    rho = sim$rho,
    `sigma_exp[1]` = 1,
    `sigma_exp[2]` = 1
  )

  want <- c(CORE, EXTRA)

  sm_raw <- fit_data$summary
  if (is.null(sm_raw) || stat != "ok") {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = stat,
      n_div = n_div, max_rhat = max_rhat
    )
    next
  }

  sm <- as.data.frame(sm_raw)
  sm$param <- rownames(sm_raw)
  rownames(sm) <- NULL
  sm <- sm[sm$param %in% want, , drop = FALSE]

  if (nrow(sm) == 0) {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = "no_params",
      n_div = n_div, max_rhat = max_rhat
    )
    next
  }

  out <- tibble(
    model = model,
    condition_id = cid,
    rep_id = rid,
    param = sm$param,
    post_mean = sm$mean,
    post_sd = sm$sd,
    l95 = sm$`2.5%`,
    u95 = sm$`97.5%`,
    truth = truth[sm$param]
  ) |>
    mutate(
      bias = ifelse(is.na(truth), NA_real_, post_mean - truth),
      rel_bias = ifelse(is.na(truth), NA_real_,
                        ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth))),
      cover95 = ifelse(is.na(truth), NA,
                       l95 <= truth & u95 >= truth),
      n_div = n_div,
      max_rhat = max_rhat,
      status = stat
    )

  rows[[i]] <- out

  if (i %% 500 == 0) message(sprintf("... processed %d/%d", i, nrow(expected)))
}

rep_tbl <- bind_rows(rows)

# Ensure proper types
if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) {
  rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
}

write_csv(rep_tbl, file.path(RES_DIR, "summary_replications_sem.csv"))

# One-row-per-fit status table (includes missing fits)
fit_status <- rep_tbl |>
  distinct(model, condition_id, rep_id, status, n_div, max_rhat)

write_csv(fit_status, file.path(RES_DIR, "fit_status_by_replication_sem.csv"))

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

write_csv(fit_status_cond, file.path(RES_DIR, "fit_status_by_condition_sem.csv"))

message("Fit status summary:")
print(table(rep_tbl$status, useNA = "ifany"))

# Aggregate results (ok fits only, consistent with Studies 1–3)
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
    ) |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .))) |>
    left_join(fit_status_cond, by = c("condition_id", "model")) |>
    mutate(RMSE = sqrt(mean_bias^2 + coalesce(emp_sd^2, 0)))

  write_csv(cond, file.path(RES_DIR, "summary_conditions_sem.csv"))
  message("✓ wrote ", nrow(cond), " rows to results_sem/summary_conditions_sem.csv")
} else {
  message("⚠ no usable draws found in any fit.")
}
