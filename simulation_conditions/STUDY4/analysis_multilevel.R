###########################################################################
# analysis_multilevel.R — summaries for ML (minimal pooling)
# 
# FIXES APPLIED:
# 1. tau_mu truth assignment now uses explicit index matching (not positional)
# 2. mu_bar truth assignment now uses explicit index matching
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")
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
  sum(vapply(sp, function(x) as.integer(sum(x[, "divergent__"])), integer(1)))
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

# Helper to extract index from parameter names like "param[1]" or "param[2,3]"
extract_single_index <- function(param_name, pattern) {
  # pattern should be like "tau_mu" to match "tau_mu[1]", "tau_mu[2]"
  regex <- paste0("^", pattern, "\\[(\\d+)\\]$")
  m <- str_match(param_name, regex)
  if (is.na(m[1, 2])) return(NA_integer_)
  as.integer(m[1, 2])
}

# list ML fits
fits <- list.files(FITS_DIR, "^fit_(EG_MLmin|NG_MLmin)_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
if (!length(fits)) {
  message("No ML fit files in ", FITS_DIR)
  quit(status = 1)
}

sim_env <- new.env()
load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_env, inherits = FALSE)) {
    p <- file.path(DATA_DIR, sprintf("sim_dataML_cond%03d_rep%03d.rds", cid, rid))
    assign(key, safe_read(p), sim_env)
  }
  get(key, sim_env, inherits = FALSE)
}

rows <- vector("list", length(fits))
for (k in seq_along(fits)) {
  fp <- fits[k]
  m <- stringr::str_match(basename(fp), "fit_(EG_MLmin|NG_MLmin)_cond(\\d+)_rep(\\d+)\\.rds")
  model <- m[2]
  cid <- as.integer(m[3])
  rid <- as.integer(m[4])

  fit <- safe_read(fp)
  stat <- if (inherits(fit, "stanfit")) {
    if (length(fit@sim) > 0 && !is.null(fit@sim$iter) && fit@sim$iter > 0 && length(fit@sim$samples) > 0) "ok" else "empty"
  } else if (inherits(fit, "error")) "bad_rds" else "not_stanfit"

  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rows[[k]] <- tibble(model, condition_id = cid, rep_id = rid, param = NA_character_, status = "sim_missing")
    next
  }

  # truths
  Phi <- sim$phi_matrix
  rho <- sim$rho
  truth_global <- c(
    phi11 = Phi[1, 1], phi12 = Phi[1, 2], phi21 = Phi[2, 1], phi22 = Phi[2, 2],
    rho = rho
  )
  if (model == "NG_MLmin") truth_global <- c(truth_global, `sigma[1]` = 1, `sigma[2]` = 1)
  if (model == "EG_MLmin") truth_global <- c(truth_global, `sigma_exp[1]` = 1, `sigma_exp[2]` = 1)

  mu_true <- sim$true_params$mu_mat # N x 2
  tau_true <- sim$true_params$tau_mu # length 2

  n_div <- count_div(fit)
  rhat <- max_rhat(fit)

  # summarise globals
  pull_sum <- function(fit, pars) {
    if (!inherits(fit, "stanfit")) {
      return(NULL)
    }
    base <- unique(sub("\\[.*$", "", pars))
    keep <- intersect(base, fit@model_pars)
    if (!length(keep)) {
      return(NULL)
    }
    s <- try(summary(fit, pars = keep)$summary, silent = TRUE)
    if (inherits(s, "try-error")) {
      return(NULL)
    }
    df <- s |>
      as.data.frame(optional = TRUE) |>
      tibble::rownames_to_column("param")
    df |> filter(param %in% pars)
  }

  # collect global params (phi, rho, sigma/sigma_exp)
  want_glob <- c(
    "phi11", "phi12", "phi21", "phi22", "rho",
    if (model == "NG_MLmin") c("sigma[1]", "sigma[2]"),
    if (model == "EG_MLmin") c("sigma_exp[1]", "sigma_exp[2]")
  )
  sm_glob <- if (stat == "ok") pull_sum(fit, want_glob) else NULL
  tab_glob <- if (!is.null(sm_glob)) {
    sm_glob |>
      mutate(
        truth = truth_global[param],
        bias = ifelse(is.na(truth), NA_real_, mean - truth),
        rel_bias = ifelse(is.na(truth), NA_real_, ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth))),
        cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
        model = model, condition_id = cid, rep_id = rid,
        n_div = n_div, max_rhat = rhat, status = stat
      ) |>
      select(model, condition_id, rep_id, param,
        post_mean = mean, post_sd = sd, l95 = `2.5%`, u95 = `97.5%`,
        truth, bias, rel_bias, cover95, n_div, max_rhat, status
      )
  } else {
    tibble()
  }

  # unit-level mus
  tab_mu <- tibble()
  if (stat == "ok") {
    sm_mu <- pull_sum(fit, "mu")
    if (!is.null(sm_mu) && nrow(sm_mu)) {
      # parse names like "mu[3,1]" => unit=3, j=1
      idx <- stringr::str_match(sm_mu$param, "^mu\\[(\\d+),(\\d+)\\]$")
      sm_mu$unit <- as.integer(idx[, 2])
      sm_mu$j <- as.integer(idx[, 3])
      sm_mu <- sm_mu |> filter(!is.na(unit), !is.na(j))
      sm_mu$truth <- mu_true[cbind(sm_mu$unit, sm_mu$j)]
      sm_mu$parnm <- ifelse(sm_mu$j == 1, "mu[1]", "mu[2]")

      tab_mu <- sm_mu |>
        mutate(
          bias = ifelse(is.na(truth), NA_real_, mean - truth),
          rel_bias = ifelse(is.na(truth), NA_real_, ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth))),
          cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
          model = model, condition_id = cid, rep_id = rid,
          n_div = n_div, max_rhat = rhat, status = stat
        ) |>
        transmute(model, condition_id, rep_id, unit,
          param = parnm,
          post_mean = mean, post_sd = sd, l95 = `2.5%`, u95 = `97.5%`,
          truth, bias, rel_bias, cover95, n_div, max_rhat, status
        )
    }
  }

  # hyperparameters mu_bar and tau_mu (truths known)
  tab_hyp <- tibble()
  if (stat == "ok") {
    # --- mu_bar: truth = (0, 0) ---
    sm_bar <- pull_sum(fit, c("mu_bar[1]", "mu_bar[2]"))
    if (!is.null(sm_bar) && nrow(sm_bar) > 0) {
      # FIX: Use explicit index matching instead of positional assignment
      sm_bar <- sm_bar |>
        mutate(
          idx = as.integer(str_match(param, "^mu_bar\\[(\\d+)\\]$")[, 2]),
          truth = 0  # mu_bar truth is always 0 for both indices
        ) |>
        filter(!is.na(idx))
      
      tb <- sm_bar |>
        mutate(
          bias = mean - truth,
          rel_bias = ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth)),
          cover95 = (`2.5%` <= truth & `97.5%` >= truth),
          model = model, condition_id = cid, rep_id = rid, 
          n_div = n_div, max_rhat = rhat, status = stat
        ) |>
        select(model, condition_id, rep_id, param,
          post_mean = mean, post_sd = sd, l95 = `2.5%`, u95 = `97.5%`,
          truth, bias, rel_bias, cover95, n_div, max_rhat, status
        )
      tab_hyp <- bind_rows(tab_hyp, tb)
    }
    
    # --- tau_mu: truth = tau_true (condition-specific) ---
    sm_tau <- pull_sum(fit, c("tau_mu[1]", "tau_mu[2]"))
    if (!is.null(sm_tau) && nrow(sm_tau) > 0) {
      # FIX: Use explicit index matching instead of positional assignment
      sm_tau <- sm_tau |>
        mutate(
          idx = as.integer(str_match(param, "^tau_mu\\[(\\d+)\\]$")[, 2])
        ) |>
        filter(!is.na(idx)) |>
        mutate(
          truth = tau_true[idx]  # Explicit index-based lookup
        )
      
      tb <- sm_tau |>
        mutate(
          bias = mean - truth,
          rel_bias = ifelse(abs(truth) < .Machine$double.eps, bias, bias / abs(truth)),
          cover95 = (`2.5%` <= truth & `97.5%` >= truth),
          model = model, condition_id = cid, rep_id = rid, 
          n_div = n_div, max_rhat = rhat, status = stat
        ) |>
        select(model, condition_id, rep_id, param,
          post_mean = mean, post_sd = sd, l95 = `2.5%`, u95 = `97.5%`,
          truth, bias, rel_bias, cover95, n_div, max_rhat, status
        )
      tab_hyp <- bind_rows(tab_hyp, tb)
    }
  }

  rows[[k]] <- bind_rows(tab_glob, tab_mu, tab_hyp)
  if (k %% 200 == 0) message(sprintf("... processed %d/%d fits", k, length(fits)))
}

rep_tbl <- bind_rows(rows)

# Ensure unit column exists even when no unit-level summaries were produced
if (!"unit" %in% names(rep_tbl)) rep_tbl$unit <- NA_integer_

out_rep <- file.path(RES_DIR, "summary_replications_ml.csv")
write_csv(rep_tbl, out_rep)

message("Fit status summary (ML):")
print(table(Status = rep_tbl$status, Model = rep_tbl$model, useNA = "ifany"))

# Check for sigma_exp recovery (diagnostic)
sigma_exp_count <- sum(rep_tbl$param %in% c("sigma_exp[1]", "sigma_exp[2]"), na.rm = TRUE)
message(sprintf("sigma_exp parameters recovered: %d rows", sigma_exp_count))

# aggregate to condition level
good <- rep_tbl$status == "ok" & !is.na(rep_tbl$param)
cond_tbl <- rep_tbl |>
  filter(good) |>
  mutate(is_unit = !is.na(unit)) |>
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
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

out_cond <- file.path(RES_DIR, "summary_conditions_ml.csv")
write_csv(cond_tbl, out_cond)
message("✓ wrote ML summaries:\n  - ", out_rep, "\n  - ", out_cond)
