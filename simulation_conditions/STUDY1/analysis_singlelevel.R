###########################################################################
# analysis_singlelevel.R  - Single-CPU analysis of fit summaries
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ── Folders ---------------------------------------------------------------
if (!requireNamespace("this.path", quietly = TRUE)) {
  message("Package 'this.path' not found. Assuming current working directory.")
  BASE_DIR <- getwd()
} else {
  tryCatch({
    BASE_DIR <- this.path::this.dir()
  }, error = function(e) {
    message("this.path::this.dir() failed. Assuming current working directory.")
    BASE_DIR <<- getwd()
  })
}

DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, FALSE, TRUE)

# ── Helpers ---------------------------------------------------------------
safe_read <- function(p) tryCatch(readRDS(p), error = function(e) e)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Parameter sets --------------------------------------------------------
CORE <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
SIGMA_NG <- c("sigma[1]", "sigma[2]")
EXTRA_SG <- c("omega[1]", "omega[2]", "alpha[1]", "alpha[2]")

# ── Enumerate fits --------------------------------------------------------
fits <- list.files(FITS_DIR, "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$",
                   full.names = TRUE)

if (length(fits) == 0) {
  message("No fit files found in ", FITS_DIR)
  if (sys.nframe() == 0) q(status = 1) else return(invisible())
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

message("Analysing ", length(fits), " fits...")
rows <- vector("list", length(fits))

for (i in seq_along(fits)) {
  fp <- fits[i]
  m <- str_match(basename(fp), "fit_(SG|NG)_cond(\\d+)_rep(\\d+)\\.rds")
  model <- m[2]
  cid <- as.integer(m[3])
  rid <- as.integer(m[4])

  fit_data <- safe_read(fp)
  
  # Handle read errors
  if (inherits(fit_data, "error")) {
    rows[[i]] <- tibble(model, condition_id = cid, rep_id = rid,
                        param = NA_character_, status = "bad_rds")
    next
  }
  
  # Extract status from the new format
  stat <- fit_data$meta$status %||% "unknown"
  
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

  # Extract summary from fit data
  sm <- fit_data$summary
  diag <- fit_data$diagnostics
  
  if (is.null(sm) || stat != "ok") {
    rows[[i]] <- tibble(
      model, condition_id = cid, rep_id = rid,
      param = NA_character_, status = stat,
      n_div = diag$n_div %||% NA_integer_,
      max_rhat = diag$max_rhat %||% NA_real_
    )
    next
  }

  # Filter to wanted parameters
  sm <- sm |> filter(param %in% want)
  
  if (nrow(sm) == 0) {
    rows[[i]] <- tibble(model, condition_id = cid, rep_id = rid,
                        param = NA_character_, status = "no_params")
    next
  }

  rows[[i]] <- sm |>
    mutate(
      truth = truth[param],
      bias = ifelse(is.na(truth), NA_real_, mean - truth),
      rel_bias = ifelse(is.na(truth), NA_real_,
                        ifelse(abs(truth) < .Machine$double.eps,
                               bias, bias / abs(truth))),
      cover95 = ifelse(is.na(truth), NA, (`2.5%` <= truth & `97.5%` >= truth)),
      model = model, condition_id = cid, rep_id = rid,
      n_div = diag$n_div %||% NA_integer_,
      max_rhat = diag$max_rhat %||% NA_real_,
      status = stat
    ) |>
    select(model, condition_id, rep_id, param,
           post_mean = mean, post_sd = sd,
           l95 = `2.5%`, u95 = `97.5%`,
           truth, bias, rel_bias, cover95,
           n_div, max_rhat, status)

  if (i %% 500 == 0) message(sprintf("... processed %d/%d", i, length(fits)))
}

rep_tbl <- bind_rows(rows)

# Ensure proper types
if ("cover95" %in% names(rep_tbl) && is.logical(rep_tbl$cover95)) {
  rep_tbl$cover95 <- as.numeric(rep_tbl$cover95)
}

write_csv(rep_tbl, file.path(RES_DIR, "summary_replications.csv"))

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

  write_csv(cond, file.path(RES_DIR, "summary_conditions.csv"))
  message("Wrote ", nrow(cond), " condition rows to summary_conditions.csv")
} else {
  message("No usable draws found in any fit.")
}
