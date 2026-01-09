###########################################################################
# analysis_singlelevel.R – Study 2 (Gamma extension)
#
# Reads fit artifacts from /fits and simulation truth from /data, then writes
#   • results/summary_replications.csv
#   • results/summary_conditions.csv
#
# This version supports the two fitted models:
#   • NG (Normal margins + Gaussian copula)
#   • GG (Gamma margins + Gaussian copula; shape estimated)
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── folders --------------------------------------------------------------
get_base_dir <- function() {
  if (requireNamespace("this.path", quietly = TRUE)) {
    return(this.path::this.dir())
  }
  ofile <- NULL
  try(ofile <- sys.frame(1)$ofile, silent = TRUE)
  if (!is.null(ofile)) {
    return(normalizePath(dirname(ofile)))
  }
  getwd()
}

BASE_DIR <- get_base_dir()
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR  <- file.path(BASE_DIR, "results")

dir.create(RES_DIR, FALSE, TRUE)
message("Analyzing results in ", BASE_DIR)

# ── helpers --------------------------------------------------------------
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

count_div <- function(fit) {
  # Lightweight fit artifact saved by fit_models.R in "summary" mode
  if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    return(as.integer(fit$n_div %||% NA_integer_))
  }
  if (!inherits(fit, "stanfit")) return(NA_integer_)
  sp <- try(get_sampler_params(fit, inc_warmup = FALSE), silent = TRUE)
  if (inherits(sp, "try-error") || length(sp) == 0) return(NA_integer_)
  sum(vapply(sp, \(x) as.integer(sum(x[, "divergent__"])), 0L))
}

max_rhat <- function(fit) {
  if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    s <- fit$summary
    if (is.data.frame(s) && ("Rhat" %in% names(s))) {
      r <- suppressWarnings(as.numeric(s$Rhat))
      r <- r[is.finite(r)]
      return(if (length(r)) max(r) else NA_real_)
    }
    return(NA_real_)
  }
  if (!inherits(fit, "stanfit")) return(NA_real_)
  s <- try(summary(fit)$summary, silent = TRUE)
  if (inherits(s, "try-error") || !is.matrix(s)) return(NA_real_)
  r <- suppressWarnings(as.numeric(s[, "Rhat"]))
  r <- r[is.finite(r)]
  if (!length(r)) return(NA_real_)
  max(r)
}

extract_fit_summary <- function(fit, pars) {
  if (is.list(fit) && identical(fit$type, "stan_fit_summary_v1")) {
    s <- fit$summary
    if (!is.data.frame(s)) return(NULL)
    s <- s |> filter(.data$param %in% pars)
    if (!nrow(s)) return(NULL)
    return(s)
  }

  if (!inherits(fit, "stanfit")) return(NULL)
  keep <- intersect(pars, fit@model_pars)
  if (!length(keep)) return(NULL)

  sm <- rstan::summary(fit, pars = keep)$summary
  df <- as.data.frame(sm)
  df$param <- rownames(df)
  rownames(df) <- NULL

  # Normalise column names to match the lightweight summary schema
  # (post_mean, post_sd, l95, u95, n_eff, Rhat)
  df <- df |>
    transmute(
      param = .data$param,
      post_mean = .data$mean,
      post_sd   = .data$sd,
      l95       = .data$`2.5%`,
      u95       = .data$`97.5%`,
      n_eff     = .data$n_eff,
      Rhat      = .data$Rhat
    )
  df
}

# ── enumerate fit artifacts ---------------------------------------------
fit_files <- list.files(FITS_DIR, "^fit_(NG|GG)_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
if (!length(fit_files)) {
  message("No fit files found in ", FITS_DIR)
  return(invisible())
}

meta <- stringr::str_match(basename(fit_files), "^fit_(NG|GG)_cond(\\d+)_rep(\\d+)\\.rds$")
meta <- as.data.frame(meta, stringsAsFactors = FALSE)
colnames(meta) <- c("file", "model", "condition_id", "rep_id")
meta$condition_id <- as.integer(meta$condition_id)
meta$rep_id <- as.integer(meta$rep_id)
meta$path <- fit_files

# ── parameter sets -------------------------------------------------------
CORE_PARS <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
WANT <- list(
  NG = c(CORE_PARS, "sigma[1]", "sigma[2]"),
  GG = c(CORE_PARS, "sigma_gam[1]", "sigma_gam[2]", "shape_gam")
)

# ── build replication-level table ---------------------------------------
rows <- vector("list", length = nrow(meta))

for (i in seq_len(nrow(meta))) {
  m <- meta[i, ]
  fit_obj <- safe_read(m$path)

  status <- "ok"
  if (inherits(fit_obj, "error")) status <- "read_error"

  n_div <- count_div(fit_obj)
  rhat  <- max_rhat(fit_obj)

  # Load corresponding simulation truth
  sim_path <- file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", m$condition_id, m$rep_id))
  sim_obj <- safe_read(sim_path)

  if (inherits(sim_obj, "error")) {
    status <- paste0(status, "|missing_sim")
    truth <- list(mu = c(0, 0), phi = matrix(NA_real_, 2, 2), rho_eff = NA_real_, shape = NA_real_)
  } else {
    phi_true <- sim_obj$phi_matrix

    # Determine effective rho under mirroring
    mir1 <- isTRUE(sim_obj$true_params$margin1$mirror)
    mir2 <- isTRUE(sim_obj$true_params$margin2$mirror)
    s1 <- if (mir1) -1 else 1
    s2 <- if (mir2) -1 else 1
    rho_eff <- (s1 * s2) * sim_obj$rho

    # Gamma shape truth (if available)
    shape_true <- NA_real_
    if (!is.null(sim_obj$true_params$margin1$type) && sim_obj$true_params$margin1$type == "gamma") {
      shape_true <- suppressWarnings(as.numeric(sim_obj$true_params$margin1$shape))
    }

    truth <- list(mu = c(0, 0), phi = phi_true, rho_eff = rho_eff, shape = shape_true)
  }

  # Truth map for each parameter
  truth_map <- function(param) {
    if (param == "mu[1]") return(truth$mu[1])
    if (param == "mu[2]") return(truth$mu[2])
    if (param == "phi11") return(truth$phi[1, 1])
    if (param == "phi12") return(truth$phi[1, 2])
    if (param == "phi21") return(truth$phi[2, 1])
    if (param == "phi22") return(truth$phi[2, 2])
    if (param == "rho")   return(truth$rho_eff)

    if (param %in% c("sigma[1]", "sigma[2]")) return(1)
    if (param %in% c("sigma_gam[1]", "sigma_gam[2]")) return(1)
    if (param == "shape_gam") return(truth$shape)

    NA_real_
  }

  summ <- extract_fit_summary(fit_obj, WANT[[m$model]])
  if (is.null(summ)) {
    rows[[i]] <- tibble(
      condition_id = m$condition_id,
      rep_id       = m$rep_id,
      model        = m$model,
      status       = status,
      n_div        = n_div,
      max_rhat     = rhat,
      param        = NA_character_,
      post_mean    = NA_real_,
      post_sd      = NA_real_,
      l95          = NA_real_,
      u95          = NA_real_,
      n_eff        = NA_real_,
      truth        = NA_real_
    )
    next
  }

  summ <- summ |>
    mutate(
      truth = vapply(.data$param, truth_map, numeric(1))
    )

  rows[[i]] <- tibble(
    condition_id = m$condition_id,
    rep_id       = m$rep_id,
    model        = m$model,
    status       = status,
    n_div        = n_div,
    max_rhat     = rhat
  ) |>
    crossing(summ)
}

rep_df <- bind_rows(rows) |>
  mutate(
    # Bias / relative bias
    bias = .data$post_mean - .data$truth,
    rel_bias = dplyr::if_else(is.na(.data$truth), NA_real_,
      dplyr::if_else(abs(.data$truth) < 1e-8, .data$bias, .data$bias / abs(.data$truth))
    ),
    cover95 = dplyr::if_else(is.na(.data$truth), NA, (.data$l95 <= .data$truth) & (.data$u95 >= .data$truth))
  )

# Drop placeholder rows (failed extraction) from downstream condition summaries
rep_df_valid <- rep_df |> filter(!is.na(.data$param))

if (!nrow(rep_df_valid)) {
  message("⚠ no usable draws found in any fit.")
  readr::write_csv(rep_df, file.path(RES_DIR, "summary_replications.csv"))
  return(invisible())
}

# ── condition-level summaries ------------------------------------------
cond_df <- rep_df_valid |>
  group_by(condition_id, model, param) |>
  summarise(
    n_reps = n(),
    truth  = dplyr::first(truth),
    mean_post_mean = mean(post_mean, na.rm = TRUE),
    mean_post_sd   = mean(post_sd, na.rm = TRUE),
    mean_bias      = mean(bias, na.rm = TRUE),
    mean_rel_bias  = mean(rel_bias, na.rm = TRUE),
    coverage_95    = mean(as.numeric(cover95), na.rm = TRUE),
    emp_sd         = sd(post_mean, na.rm = TRUE),
    mean_n_div     = mean(n_div, na.rm = TRUE),
    mean_max_rhat  = mean(max_rhat, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    sd_bias = mean_post_sd - dplyr::if_else(is.na(emp_sd), 0, emp_sd)
  )

# ── write outputs --------------------------------------------------------
readr::write_csv(rep_df,  file.path(RES_DIR, "summary_replications.csv"))
readr::write_csv(cond_df, file.path(RES_DIR, "summary_conditions.csv"))

# Small console summary
fit_status <- rep_df |>
  distinct(condition_id, rep_id, model, status) |>
  count(model, status) |>
  arrange(model, status)
message("Fit status summary:")
print(fit_status)
