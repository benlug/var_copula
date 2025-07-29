#!/usr/bin/env Rscript
###########################################################################
# analysis_singlelevel.R  – sequential, robust summariser (fixed type)
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr);  library(tidyr)
  library(readr);  library(stringr)
})

# ── folders --------------------------------------------------------------
BASE_DIR <- this.path::this.dir()
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR  <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, FALSE, TRUE)

# ── helpers --------------------------------------------------------------
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%`    <- function(a, b) if (!is.null(a)) a else b

count_div <- function(fit) {
  if (!inherits(fit, "stanfit")) return(NA_integer_)
  tot <- sum(vapply(rstan::get_sampler_params(fit, FALSE),
                    \(x) sum(x[, "divergent__"]), numeric(1)))
  as.integer(tot)
}

quick_summary <- function(fit, want_full) {
  base <- unique(sub("\\[.*$", "", want_full))
  keep <- intersect(base, fit@model_pars)
  if (!length(keep)) return(NULL)

  rstan::summary(fit, pars = keep)$summary |>
    as.data.frame(optional = TRUE) |>
    tibble::rownames_to_column("param") |>
    dplyr::filter(param %in% want_full)
}

# ── parameter sets -------------------------------------------------------
CORE      <- c("mu[1]","mu[2]","phi11","phi12","phi21","phi22","rho")
SIGMA_NG  <- c("sigma[1]","sigma[2]")
EXTRA_SG  <- c("omega[1]","omega[2]","delta[1]","delta[2]")

# ── fit files ------------------------------------------------------------
fit_files <- list.files(FITS_DIR,
                        "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$",
                        full.names = TRUE)
if (!length(fit_files))
  stop("No fit_*.rds files found in ", FITS_DIR)

# ── simulation cache -----------------------------------------------------
sim_cache <- new.env(parent = emptyenv())
load_sim <- function(cid, rid) {
  key <- sprintf("%03d_%03d", cid, rid)
  if (!exists(key, sim_cache, inherits = FALSE)) {
    p <- file.path(DATA_DIR,
                   sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid))
    assign(key, safe_read(p), sim_cache)
  }
  get(key, sim_cache, inherits = FALSE)
}

# ── sequential loop ------------------------------------------------------
message("⏳ analysing ", length(fit_files), " fits on ONE core ...")
rep_rows <- vector("list", length(fit_files))

for (k in seq_along(fit_files)) {

  fp <- fit_files[k]
  m  <- stringr::str_match(basename(fp),
                           "fit_(SG|NG)_cond(\\d+)_rep(\\d+)\\.rds")
  model <- m[2];  cid <- as.integer(m[3]);  rid <- as.integer(m[4])

  fit_obj <- safe_read(fp)
  status <- if (inherits(fit_obj, "stanfit")) {
    it <- tryCatch(fit_obj@sim$iter, error = \(e) NA_integer_)
    if (is.na(it))      "not_iter"
    else if (it == 0L)  "empty"
    else                "ok"
  } else if (inherits(fit_obj, "error")) "bad_rds" else "not_stanfit"

  sim <- load_sim(cid, rid)
  if (inherits(sim, "error")) {
    rep_rows[[k]] <- tibble(model, condition_id = cid, rep_id = rid,
                            param = NA_character_, status,
                            n_div = NA_integer_, note = "sim missing")
    next
  }

  # --- ground truth ------------------------------------------------------
  thr <- with(sim, c(
    `mu[1]`=0, `mu[2]`=0,
    phi11=phi_matrix[1,1], phi12=phi_matrix[1,2],
    phi21=phi_matrix[2,1], phi22=phi_matrix[2,2],
    rho=rho))
  if (model=="NG")
    thr <- c(thr, `sigma[1]`=1, `sigma[2]`=1)
  if (model=="SG") {
    a1 <- sim$true_params$margin1$alpha %||% 0
    a2 <- sim$true_params$margin2$alpha %||% 0
    d1 <- a1 / sqrt(1 + a1^2);  d2 <- a2 / sqrt(1 + a2^2)
    thr <- c(thr,
             `omega[1]`=1, `omega[2]`=1,
             `delta[1]`=d1, `delta[2]`=d2)
  }

  want <- c(CORE,
            if (model=="NG") SIGMA_NG,
            if (model=="SG") EXTRA_SG)

  sm <- NULL
  if (status=="ok") sm <- quick_summary(fit_obj, want)

  if (is.null(sm)) {
    rep_rows[[k]] <- tibble(model, condition_id = cid, rep_id = rid,
                            param = NA_character_, status,
                            n_div = count_div(fit_obj))
  } else {
    rep_rows[[k]] <- sm |>
      mutate(
        truth     = thr[param],
        bias      = mean - truth,
        rel_bias  = ifelse(abs(truth) < .Machine$double.eps,
                           bias, bias / abs(truth)),
        cover95   = (`2.5%` <= truth & `97.5%` >= truth),
        model     = model,
        condition_id = cid,
        rep_id    = rid,
        n_div     = count_div(fit_obj),
        status    = status
      ) |>
      select(model, condition_id, rep_id, param,
             post_mean = mean, post_sd = sd,
             l95 = `2.5%`, u95 = `97.5%`,
             truth, bias, rel_bias, cover95, n_div, status)
  }
}

rep_tbl <- dplyr::bind_rows(rep_rows)
readr::write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))
print(table(rep_tbl$status, useNA = "ifany"))

# --- aggregate -----------------------------------------------------------
usable <- with(rep_tbl, status=="ok" & !is.na(param))
if (any(usable)) {
  agg <- rep_tbl |>
    dplyr::filter(usable) |>
    dplyr::group_by(condition_id, model, param) |>
    dplyr::summarise(
      mean_rel_bias = mean(rel_bias),
      coverage_95   = mean(cover95),
      mean_post_sd  = mean(post_sd),
      emp_sd        = sd(post_mean),
      sd_bias       = mean_post_sd - emp_sd,
      mean_n_div    = mean(n_div),
      .groups = "drop"
    )
  readr::write_csv(agg, file.path(RES_DIR, "summary_conditions.csv"))
  message("✓ wrote ", nrow(agg), " condition rows")
} else {
  message("⚠ no usable draws – all fits empty / corrupt / missing")
}
