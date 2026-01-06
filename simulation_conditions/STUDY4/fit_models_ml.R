###########################################################################
# fit_models_ml.R — Fit multilevel VAR(1) copula models (NG_MLmin, EG_MLmin)
#
# Expects simulation files:
#   data/sim_dataML_cond%03d_rep%03d.rds
#
# Writes fit objects:
#   fits/fit_{EG_MLmin|NG_MLmin}_cond%03d_rep%03d.rds
#
# Resume / design safety:
#   - If data/sim_conditions_ml.rds exists, only those condition_ids are processed.
#   - Existing fit files are skipped (unless overwrite=TRUE).
#   - Deterministic Stan seeds per (condition_id, rep_id, model).
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Fit VAR(1) copula models for multilevel (minimal pooling) simulation outputs
#'
#' @param data_dir Directory containing sim_dataML_*.rds files
#' @param fits_dir Directory to write fit_*.rds files
#' @param stan_dir Directory containing Stan model files
#' @param results_dir (Optional) Directory for any auxiliary outputs (not required)
#' @param chains Number of chains
#' @param iter Total iterations per chain
#' @param warmup Warmup iterations per chain
#' @param adapt_delta Stan adapt_delta
#' @param max_treedepth Stan max_treedepth
#' @param cores_outer Intended outer parallelism level (used only on unix via mclapply)
#' @param start_condition Resume: minimum condition_id to process
#' @param start_rep Resume: minimum rep_id for start_condition
#' @param overwrite If TRUE, refit even when fit files exist
#'
#' @return Invisibly, a data frame of tasks processed.
fit_var1_copula_models_ml <- function(
    data_dir,
    fits_dir,
    stan_dir,
    results_dir = NULL,
    chains = 4,
    iter = 3000,
    warmup = 1500,
    adapt_delta = 0.90,
    max_treedepth = 10,
    cores_outer = 1,
    start_condition = 1,
    start_rep = 1,
    overwrite = FALSE
) {
  dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(results_dir)) dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  # Keep chain parallelism OFF when running outer parallelism to avoid oversubscription.
  # (If you prefer parallel chains, set options(mc.cores = chains) and keep cores_outer = 1.)
  options(mc.cores = 1)
  rstan_options(auto_write = TRUE)

  eg_file <- file.path(stan_dir, "model_EG_ml_min.stan")
  ng_file <- file.path(stan_dir, "model_NG_ml_min.stan")
  if (!file.exists(eg_file)) stop("Missing Stan file: ", eg_file)
  if (!file.exists(ng_file)) stop("Missing Stan file: ", ng_file)

  message("Compiling Stan models …")
  sm_EG <- rstan::stan_model(eg_file)
  sm_NG <- rstan::stan_model(ng_file)

  sim_files <- list.files(data_dir, "^sim_dataML_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(sim_files)) {
    stop("No simulation files found in ", data_dir)
  }

  # Restrict to active design conditions (if design file exists)
  allowed_cids <- NULL
  design_file <- file.path(data_dir, "sim_conditions_ml.rds")
  if (file.exists(design_file)) {
    design <- try(readRDS(design_file), silent = TRUE)
    if (!inherits(design, "try-error") && is.data.frame(design) && "condition_id" %in% names(design)) {
      allowed_cids <- sort(unique(as.integer(design$condition_id)))
      allowed_cids <- allowed_cids[is.finite(allowed_cids)]
    }
  }

  tasks <- tibble::tibble(file = sim_files) %>%     mutate(
      base = basename(file),
      condition_id = as.integer(str_match(base, "sim_dataML_cond(\\d+)_rep(\\d+)\\.rds")[, 2]),
      rep_id       = as.integer(str_match(base, "sim_dataML_cond(\\d+)_rep(\\d+)\\.rds")[, 3])
    ) %>%     filter(!is.na(condition_id), !is.na(rep_id))

  if (!is.null(allowed_cids)) {
    extra <- setdiff(sort(unique(tasks$condition_id)), allowed_cids)
    if (length(extra)) {
      message(
        "Design file found (", design_file, ").\n",
        "Ignoring simulation files for condition_id(s) not in the design: ",
        paste(extra, collapse = ", ")
      )
    }
    tasks <- tasks %>% filter(condition_id %in% allowed_cids)
  }

  # Resume filters
  tasks <- tasks %>%     filter(condition_id > start_condition | (condition_id == start_condition & rep_id >= start_rep)) %>%     arrange(condition_id, rep_id)

  # Pre-compute output paths and drop tasks already finished (unless overwrite=TRUE)
  tasks <- tasks %>%     mutate(
      out_ng  = file.path(fits_dir, sprintf("fit_NG_MLmin_cond%03d_rep%03d.rds", condition_id, rep_id)),
      out_eg  = file.path(fits_dir, sprintf("fit_EG_MLmin_cond%03d_rep%03d.rds", condition_id, rep_id)),
      need_ng = overwrite | !file.exists(out_ng),
      need_eg = overwrite | !file.exists(out_eg),
      need_any = need_ng | need_eg
    ) %>%     filter(need_any)

  if (!nrow(tasks)) {
    message("No tasks to run after applying filters (resume + already-fitted).")
    return(invisible(tasks))
  }

  fit_one <- function(task_row) {
    cid <- task_row$condition_id
    rid <- task_row$rep_id

    sim <- readRDS(task_row$file)

    # Build y as array[N] matrix[T,2] for Stan
    panel <- sim$data
    N <- sim$N
    Tn <- sim$T
    y_list <- vector("list", N)
    for (i in seq_len(N)) {
      df_i <- panel[panel$i == i, , drop = FALSE]
      df_i <- df_i[order(df_i$t), , drop = FALSE]
      y_list[[i]] <- as.matrix(df_i[, c("y1", "y2")])
    }

    # EG skew direction (shared across units)
    skew_dir <- sim$true_params$skew_dir %||% c(1, 1)
    skew_dir <- as.numeric(skew_dir[1:2])
    skew_dir <- ifelse(skew_dir < 0, -1, 1)

    stan_data_ng <- list(N = N, T = Tn, y = y_list)
    stan_data_eg <- list(N = N, T = Tn, y = y_list, skew_direction = as.vector(skew_dir))

    # Deterministic seeds per (condition, replication, model)
    base_seed <- 910000L + as.integer(cid) * 10000L + as.integer(rid)

    # --- NG fit ---
    if (isTRUE(task_row$need_ng)) {
      fit_ng <- try(
        rstan::sampling(
          object = sm_NG,
          data = stan_data_ng,
          chains = chains,
          iter = iter,
          warmup = warmup,
          seed = base_seed + 1L,
          refresh = 0,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
        ),
        silent = TRUE
      )
      if (inherits(fit_ng, "try-error")) {
        saveRDS(structure(list(message = as.character(fit_ng)), class = "error"), task_row$out_ng)
      } else {
        saveRDS(fit_ng, task_row$out_ng)
      }
    }

    # --- EG fit ---
    if (isTRUE(task_row$need_eg)) {
      fit_eg <- try(
        rstan::sampling(
          object = sm_EG,
          data = stan_data_eg,
          chains = chains,
          iter = iter,
          warmup = warmup,
          seed = base_seed + 2L,
          refresh = 0,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
        ),
        silent = TRUE
      )
      if (inherits(fit_eg, "try-error")) {
        saveRDS(structure(list(message = as.character(fit_eg)), class = "error"), task_row$out_eg)
      } else {
        saveRDS(fit_eg, task_row$out_eg)
      }
    }

    tibble::tibble(condition_id = cid, rep_id = rid, did_ng = isTRUE(task_row$need_ng), did_eg = isTRUE(task_row$need_eg))
  }

  # Outer parallelism (best-effort): only use forked mclapply on unix.
  if (cores_outer > 1 && .Platform$OS.type == "unix") {
    message("Running outer parallel fits with mclapply (cores_outer = ", cores_outer, ") …")
    message("(Progress bar disabled in parallel mode.)")
    res <- parallel::mclapply(
      seq_len(nrow(tasks)),
      function(i) fit_one(tasks[i, , drop = FALSE]),
      mc.cores = cores_outer
    )
    res <- dplyr::bind_rows(res)
  } else {
    if (cores_outer > 1) {
      message("cores_outer > 1 requested, but OS is not unix; running sequentially.")
    }
    message("Fitting ", nrow(tasks), " dataset(s) …")
    pb <- utils::txtProgressBar(min = 0, max = nrow(tasks), style = 3)
    res_list <- vector("list", nrow(tasks))
    for (k in seq_len(nrow(tasks))) {
      res_list[[k]] <- fit_one(tasks[k, , drop = FALSE])
      utils::setTxtProgressBar(pb, k)
    }
    close(pb)
    cat("\n")
    res <- dplyr::bind_rows(res_list)
  }

  message("✓ Finished fitting ML models.")
  invisible(res)
}
