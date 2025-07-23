###########################################################################
# investigate_failed_fits.R  –  bullet‑proof version
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(rstan)
})

## -- folders -------------------------------------------------------------
BASE_DIR <- this.path::this.dir()
DATA_DIR <- file.path(BASE_DIR, "../data")
FITS_DIR <- file.path(BASE_DIR, "../fits")
RESULTS_DIR <- file.path(BASE_DIR, "../results")
STAN_DIR <- file.path(BASE_DIR, "../stan")
dir.create(RESULTS_DIR, FALSE, TRUE)

## ------------------------------------------------------------------------
## helper: classify one fit file → 1‑row tibble ---------------------------
## ------------------------------------------------------------------------
inspect_fit_file <- function(path, model, cid, rid) {
  safe_row <- function(stat = "missing", it = NA_integer_, div = NA_integer_) {
    tibble(
      model = model, condition_id = cid, rep_id = rid,
      status = stat, iter = it, n_div = div
    )
  }

  if (!file.exists(path)) {
    return(safe_row("missing"))
  }

  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) {
    return(safe_row("corrupt"))
  }

  if (!inherits(obj, "stanfit")) {
    return(safe_row("corrupt"))
  }

  it <- tryCatch(obj@sim$iter, error = function(e) NA_integer_)
  if (length(it) != 1 || is.na(it) || it == 0) {
    return(safe_row("empty", 0))
  }

  n_div <- tryCatch(
    sum(vapply(
      rstan::get_sampler_params(obj, FALSE),
      function(x) sum(x[, "divergent__"]), 0
    )),
    error = function(e) NA_integer_
  )
  safe_row("ok", it, n_div)
}

## ------------------------------------------------------------------------
## 1. expected SG / NG grid ----------------------------------------------
## ------------------------------------------------------------------------
data_files <- list.files(DATA_DIR,
  "^sim_data_cond(\\d+)_rep(\\d+)\\.rds$",
  full.names = TRUE
)
stopifnot(length(data_files) > 0)

data_df <- tibble(
  condition_id = as.integer(str_match(basename(data_files), "cond(\\d+)")[, 2]),
  rep_id       = as.integer(str_match(basename(data_files), "rep(\\d+)")[, 2])
)

exp_grid <- tidyr::crossing(data_df, model = c("SG", "NG")) |>
  mutate(
    fit_path = file.path(
      FITS_DIR,
      sprintf("fit_%s_cond%03d_rep%03d.rds", model, condition_id, rep_id)
    )
  )

## ------------------------------------------------------------------------
## 2. iterate with progress bar -------------------------------------------
## ------------------------------------------------------------------------
pb <- txtProgressBar(0, nrow(exp_grid), style = 3)
res_rows <- vector("list", nrow(exp_grid))

for (i in seq_len(nrow(exp_grid))) {
  row <- exp_grid[i, ]
  res_rows[[i]] <- tryCatch(
    inspect_fit_file(
      row$fit_path, row$model,
      row$condition_id, row$rep_id
    ),
    error = function(e) {
      tibble(
        model = row$model,
        condition_id = row$condition_id,
        rep_id = row$rep_id,
        status = "inspect_error",
        iter = NA_integer_, n_div = NA_integer_
      )
    }
  )
  setTxtProgressBar(pb, i)
}
close(pb)

fit_status_tbl <- bind_rows(res_rows)

## ------------------------------------------------------------------------
## 3. optional: merge log lines / ERROR_*.rds -----------------------------
## ------------------------------------------------------------------------
log_tbl <- {
  lp <- file.path(FITS_DIR, "stan_fit.log")
  if (file.exists(lp)) {
    read_lines(lp) |>
      as_tibble() |>
      filter(str_detect(value, "FAIL")) |>
      mutate(
        model = str_match(value, "(SG|NG)")[, 2],
        condition_id = as.integer(str_match(value, "cond(\\d+)")[, 2]),
        rep_id = as.integer(str_match(value, "rep(\\d+)")[, 2]),
        log_msg = value
      ) |>
      select(model, condition_id, rep_id, log_msg)
  } else {
    tibble(
      model = character(), condition_id = integer(),
      rep_id = integer(), log_msg = character()
    )
  }
}

err_tbl <- {
  ef <- list.files(FITS_DIR,
    "^ERROR_(SG|NG)_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(ef)) {
    tibble(
      model = str_match(basename(ef), "ERROR_(SG|NG)")[, 2],
      condition_id = as.integer(str_match(basename(ef), "cond(\\d+)")[, 2]),
      rep_id = as.integer(str_match(basename(ef), "rep(\\d+)")[, 2]),
      has_error_rds = TRUE,
      err_file = ef
    )
  } else {
    tibble(
      model = character(), condition_id = integer(),
      rep_id = integer(), has_error_rds = logical(), err_file = character()
    )
  }
}

## ------------------------------------------------------------------------
## 4. final failure table -------------------------------------------------
## ------------------------------------------------------------------------
fail_tbl <- fit_status_tbl |>
  left_join(log_tbl, by = c("model", "condition_id", "rep_id")) |>
  left_join(err_tbl, by = c("model", "condition_id", "rep_id")) |>
  mutate(
    log_msg       = coalesce(log_msg, "(no log entry)"),
    has_error_rds = coalesce(has_error_rds, FALSE)
  ) |>
  filter(status != "ok") |>
  arrange(condition_id, rep_id, model)

write_csv(fail_tbl, file.path(RESULTS_DIR, "failed_fits.csv"))

message(
  "✓  ", nrow(fail_tbl),
  " non‑OK fits → results/failed_fits.csv"
)

## ------------------------------------------------------------------------
## 5. helper to re‑run one dataset with verbose Stan output ---------------
## ------------------------------------------------------------------------
rerun_one_fit <- function(condition_id, rep_id,
                          model = c("SG", "NG"),
                          iter = 4000, warmup = 2000,
                          adapt_delta = 0.95,
                          max_treedepth = 12) {
  model <- match.arg(model)
  dpath <- file.path(
    DATA_DIR,
    sprintf(
      "sim_data_cond%03d_rep%03d.rds",
      condition_id, rep_id
    )
  )
  if (!file.exists(dpath)) {
    stop("Data file not found: ", dpath)
  }

  ds <- readRDS(dpath)
  stan_data <- list(
    T = ds$T,
    y = lapply(
      seq_len(ds$T),
      function(i) c(ds$data$y1[i], ds$data$y2[i])
    )
  )

  sfile <- if (model == "SG") {
    file.path(STAN_DIR, "model_SNG_sl.stan")
  } else {
    file.path(STAN_DIR, "model_NG_sl.stan")
  }
  mod <- rstan::stan_model(sfile)

  sampling(mod,
    data = stan_data,
    chains = 4,
    iter = iter,
    warmup = warmup,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = 1000, verbose = TRUE
  )
}

View(fail_tbl)
# rerun_one_fit(condition_id = 1, rep_id = 8, model = "SG")


# fit <- readRDS("simulation_conditions/12_06_25/fits/fit_SG_cond001_rep008.rds")
