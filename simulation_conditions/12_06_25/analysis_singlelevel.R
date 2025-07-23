###########################################################################
# analysis_singlelevel.R  –  Robust diagnostics & bias summary
#   • Handles fit_SG_* and fit_NG_* files
#   • Never crashes on length‑0 / NULL objects
#   • 2025‑07‑update:  adds  sigma[1:2]  (NG scale) to the core set
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ----------------- folders ----------------------------------------------
BASE_DIR <- this.path::this.dir()
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, FALSE, TRUE)

# ----------------- helpers ----------------------------------------------
safe_read <- function(path) tryCatch(readRDS(path), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

count_div <- function(fit) {
  if (!inherits(fit, "stanfit")) {
    return(NA_integer_)
  }
  sum(vapply(
    rstan::get_sampler_params(fit, FALSE),
    \(x) sum(x[, "divergent__"]), 0
  ))
}

safe_summary <- function(fit, pars) {
  if (!inherits(fit, "stanfit")) {
    return(NULL)
  }
  it <- tryCatch(fit@sim$iter, error = \(e) NA_integer_)
  if (length(it) == 0 || is.na(it) || it == 0) {
    return(NULL)
  }
  as.data.frame(summary(fit, pars = pars)$summary, optional = TRUE) |>
    tibble::rownames_to_column("param")
}

# ----------------- parameter sets ---------------------------------------
CORE <- c(
  "mu[1]", "mu[2]",
  "phi11", "phi12", "phi21", "phi22",
  "sigma[1]", "sigma[2]", # ‹‑‑ NEW (normal scale)
  "rho"
)
EXTRA <- c("omega[1]", "omega[2]", "alpha[1]", "alpha[2]") # SG only

# ----------------- gather fit files -------------------------------------
fit_files <- list.files(FITS_DIR,
  "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)
if (!length(fit_files)) {
  stop("No fit_*.rds files found in ", FITS_DIR)
}

rep_rows <- vector("list", length(fit_files))

for (k in seq_along(fit_files)) {
  fp <- fit_files[k]
  meta <- stringr::str_match(
    basename(fp),
    "fit_(SG|NG)_cond(\\d+)_rep(\\d+)\\.rds"
  )
  model <- meta[2] # "SG" or "NG"
  cid <- as.integer(meta[3])
  rid <- as.integer(meta[4])

  fit_obj <- safe_read(fp)

  ## —— status code ------------------------------------------------------
  status <- "bad_rds"
  if (inherits(fit_obj, "stanfit")) {
    it <- tryCatch(fit_obj@sim$iter, error = \(e) NA_integer_)
    status <- if (length(it) == 0 || is.na(it)) {
      "not_iter"
    } else if (it == 0) "empty" else "ok"
  } else if (!inherits(fit_obj, "error")) {
    status <- "not_stanfit"
  }

  ## —— load simulation truth -------------------------------------------
  sim_path <- file.path(
    DATA_DIR,
    sprintf("sim_data_cond%03d_rep%03d.rds", cid, rid)
  )
  if (!file.exists(sim_path)) {
    rep_rows[[k]] <- tibble(model,
      condition_id = cid, rep_id = rid,
      param = NA_character_, status,
      note = "sim missing"
    )
    next
  }
  sim <- readRDS(sim_path)

  truth <- c(
    `mu[1]` = 0, `mu[2]` = 0,
    phi11 = sim$phi_matrix[1, 1],
    phi12 = sim$phi_matrix[1, 2],
    phi21 = sim$phi_matrix[2, 1],
    phi22 = sim$phi_matrix[2, 2],
    `sigma[1]` = 1, # fixed‑‑unit residual SD under NG
    `sigma[2]` = 1,
    rho = sim$rho
  )
  if (model == "SG") {
    truth <- c(truth,
      `omega[1]` = 1, `omega[2]` = 1,
      `alpha[1]` = sim$true_params$margin1$alpha %||% 0,
      `alpha[2]` = sim$true_params$margin2$alpha %||% 0
    )
  }

  pars <- c(CORE, if (model == "SG") EXTRA)

  sm <- safe_summary(fit_obj, pars)

  if (is.null(sm)) { # «‑– no usable draws
    rep_rows[[k]] <- tibble(model,
      condition_id = cid, rep_id = rid,
      param = NA_character_, status,
      n_div = count_div(fit_obj)
    )
  } else {
    rep_rows[[k]] <- sm |>
      mutate(
        truth = truth[param],
        bias = mean - truth,
        rel_bias = ifelse(abs(truth) < .Machine$double.eps,
          bias, bias / abs(truth)
        ),
        cover95 = (`2.5%` <= truth & `97.5%` >= truth) * 1,
        model = model,
        condition_id = cid,
        rep_id = rid,
        n_div = count_div(fit_obj),
        status = status
      ) |>
      select(model, condition_id, rep_id, param,
        post_mean = mean, post_sd = sd,
        l95 = `2.5%`, u95 = `97.5%`,
        truth, bias, rel_bias, cover95,
        n_div, status
      )
  }
}

rep_tbl <- bind_rows(rep_rows)
write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))

usable <- rep_tbl$status == "ok" & !is.na(rep_tbl$param)
if (any(usable)) {
  agg <- rep_tbl |>
    filter(usable) |>
    group_by(condition_id, model, param) |>
    summarise(
      mean_rel_bias = mean(rel_bias),
      coverage_95 = mean(cover95),
      mean_post_sd = mean(post_sd),
      emp_sd = sd(post_mean),
      sd_bias = mean_post_sd - emp_sd,
      mean_n_div = mean(n_div),
      .groups = "drop"
    )
  write_csv(agg, file.path(RES_DIR, "summary_conditions.csv"))
  message(
    "✓ analysis_singlelevel.R – ", nrow(agg),
    " aggregated condition rows written."
  )
} else {
  message("⚠ No usable parameter draws found.")
}
