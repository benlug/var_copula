###########################################################################
# analysis_sem.R — Aggregation & metrics for SEM Study A/B (EI vs EL)
# Mirrors Study 2 analysis (bias, coverage, SD-bias, mcmc diagnostics). :contentReference[oaicite:3]{index=3}
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits_sem")
RES_DIR <- file.path(BASE_DIR, "results_sem")
dir.create(RES_DIR, FALSE, TRUE)

safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

count_div <- function(fit) {
  if (!inherits(fit, "stanfit")) {
    return(NA_integer_)
  }
  sp <- try(get_sampler_params(fit, FALSE), silent = TRUE)
  if (inherits(sp, "try-error") || length(sp) == 0) {
    return(NA_integer_)
  }
  sum(vapply(sp, \(x) as.integer(sum(x[, "divergent__"])), 0L))
}
max_rhat <- function(fit) {
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
    dplyr::filter(param %in% want_full)
}

# parameters to collect (same names in both Stan models)
CORE <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
EXTRA <- c("sigma_exp[1]", "sigma_exp[2]") # active-layer scales

fits <- list.files(FITS_DIR, "^fit_(EI|EL)_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
if (length(fits) == 0) {
  message("No SEM fits found.")
  if (sys.nframe() == 0) q(status = 1) else quit <- TRUE
}

# simulation cache
sim_env <- new.env()
load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_env, inherits = FALSE)) {
    p <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid))
    assign(key, safe_read(p), sim_env)
  }
  get(key, sim_env, inherits = FALSE)
}

message("⏳ analysing ", length(fits), " SEM fits on one core …")
rows <- vector("list", length(fits))

for (i in seq_along(fits)) {
  fp <- fits[i]
  m <- str_match(basename(fp), "fit_(EI|EL)_cond(\\d+)_rep(\\d+)\\.rds")
  model_code <- m[2]
  cid <- as.integer(m[3])
  rid <- as.integer(m[4])

  fit <- safe_read(fp)
  stat <- if (inherits(fit, "stanfit")) {
    if (length(fit@sim) > 0 && length(fit@sim$samples) > 0) "ok" else "empty"
  } else if (inherits(fit, "error")) "bad_rds" else "not_stanfit"

  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rows[[i]] <- tibble(
      model = model_code, condition_id = cid, rep_id = rid,
      param = NA_character_, status = "sim_missing"
    )
    next
  }

  # ground truth (standardized Exponential margins; rho at active layer)
  truth <- c(
    `mu[1]` = 0, `mu[2]` = 0,
    phi11 = sim$phi_matrix[1, 1], phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1], phi22 = sim$phi_matrix[2, 2],
    rho = sim$rho,
    `sigma_exp[1]` = 1, `sigma_exp[2]` = 1
  )

  want <- c(CORE, EXTRA)
  sm <- if (stat == "ok") quick_summary(fit, want) else NULL
  rhat <- if (stat == "ok") max_rhat(fit) else NA_real_
  ndiv <- count_div(fit)

  rows[[i]] <- if (is.null(sm)) {
    tibble(
      model = model_code, condition_id = cid, rep_id = rid,
      param = NA_character_, status = stat, n_div = ndiv, max_rhat = rhat
    )
  } else {
    sm |>
      mutate(
        truth = truth[param],
        bias = ifelse(is.na(truth), NA_real_, mean - truth),
        rel_bias = ifelse(is.na(truth), NA_real_,
          ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth))
        ),
        cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
        model = model_code, condition_id = cid, rep_id = rid,
        n_div = ndiv, max_rhat = rhat, status = stat
      ) |>
      select(model, condition_id, rep_id, param,
        post_mean = mean, post_sd = sd, l95 = `2.5%`, u95 = `97.5%`,
        truth, bias, rel_bias, cover95, n_div, max_rhat, status
      )
  }

  if (i %% 300 == 0) message(sprintf("... processed %d/%d", i, length(fits)))
}

rep_tbl <- bind_rows(rows)
if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
write_csv(rep_tbl, file.path(RES_DIR, "summary_replications_sem.csv"))

message("Fit status summary:")
print(table(Status = rep_tbl$status, Model = rep_tbl$model, useNA = "ifany"))

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
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

  write_csv(cond, file.path(RES_DIR, "summary_conditions_sem.csv"))
  message("✓ wrote ", nrow(cond), " condition rows to summary_conditions_sem.csv")
} else {
  message("⚠ no usable draws found.")
}
