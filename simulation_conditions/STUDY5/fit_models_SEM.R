#!/usr/bin/env Rscript
###########################################################################
# fit_models_SEM.R — Study 5 (SEM): Fit EI vs EL with rstan
#
# Output format (aligned with Studies 1–3):
#   Each fit is saved as a *compact* list (not a full stanfit) containing:
#     - status: "ok" | "error" | "unexpected" | ...
#     - summary: matrix from summary(fit, pars=... )$summary
#     - n_div, max_rhat
#     - model, condition_id, rep_id
#     - (optional) warnings
#
# This keeps disk usage manageable for large simulation grids and also
# supports stop/resume via file-existence checks.
#
# Parallelisation
#   Uses mclapply (forking) on Unix-like OS to avoid rstan object
#   serialisation issues (PSOCK clusters frequently break stanmodel
#   pointers). On Windows, falls back to sequential fitting.
###########################################################################

fit_sem_models <- function(data_dir,
                           fits_dir,
                           stan_dir,
                           results_dir,
                           chains = 4,
                           iter = 4000,
                           warmup = 2000,
                           adapt_delta = 0.95,
                           max_treedepth = 12,
                           cores_outer = 2,
                           start_condition = 1,
                           start_rep = 1,
                           overwrite_fits = FALSE) {

  suppressPackageStartupMessages({
    library(rstan)
    library(stringr)
    library(parallel)
    library(dplyr)
  })

  rstan::rstan_options(auto_write = TRUE)

  SEED_BASE_R    <- 7e6
  SEED_BASE_STAN <- 9e6

  side_dir <- file.path(results_dir, "init_sidecar_sem")
  dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(side_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(fits_dir, "stan_fit_sem.log")
  if (!file.exists(log_file)) try(file.create(log_file), silent = TRUE)

  # ---- compile models in main process (shared by forked workers) --------
  message("Compiling Stan models (Study 5 SEM)…")

  ei_path <- file.path(stan_dir, "sem_A_indicator_EG.stan")
  el_path <- file.path(stan_dir, "sem_B_latent_EG.stan")
  if (!file.exists(ei_path)) stop("EI model not found: ", ei_path)
  if (!file.exists(el_path)) stop("EL model not found: ", el_path)

  model_EI <- rstan::stan_model(ei_path, model_name = "SEM_Indicator_EG")
  message("  Compiled: EI")
  model_EL <- rstan::stan_model(el_path, model_name = "SEM_Latent_EG")
  message("  Compiled: EL")

  models <- list(EI = model_EI, EL = model_EL)

  # Parameters to monitor (keep output small)
  pars_keep <- c("mu", "phi11", "phi12", "phi21", "phi22", "rho", "sigma_exp")

  # ---- init helpers -----------------------------------------------------
  make_init_EI <- function(y, chain_seed) {
    set.seed(chain_seed)
    ylag <- rbind(c(0, 0), y[-nrow(y), , drop = FALSE])
    f1 <- try(lm(y[, 1] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
    f2 <- try(lm(y[, 2] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
    phi11 <- if (!inherits(f1, "try-error")) coef(f1)[2] else 0.55
    phi12 <- if (!inherits(f1, "try-error")) coef(f1)[3] else 0.10
    phi21 <- if (!inherits(f2, "try-error")) coef(f2)[2] else 0.10
    phi22 <- if (!inherits(f2, "try-error")) coef(f2)[3] else 0.25

    list(
      mu = c(0, 0),
      phi11 = max(min(phi11, 0.9), -0.9),
      phi12 = max(min(phi12, 0.9), -0.9),
      phi21 = max(min(phi21, 0.9), -0.9),
      phi22 = max(min(phi22, 0.9), -0.9),
      eta = rep(log(0.2), 2),
      rho_raw = rnorm(1, 0, 0.5),
      z_raw = matrix(0, nrow(y), 2)
    )
  }

  make_init_EL <- function(y, chain_seed) {
    set.seed(chain_seed)
    ylag <- rbind(c(0, 0), y[-nrow(y), , drop = FALSE])
    f1 <- try(lm(y[, 1] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
    f2 <- try(lm(y[, 2] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
    phi11 <- if (!inherits(f1, "try-error")) coef(f1)[2] else 0.55
    phi12 <- if (!inherits(f1, "try-error")) coef(f1)[3] else 0.10
    phi21 <- if (!inherits(f2, "try-error")) coef(f2)[2] else 0.10
    phi22 <- if (!inherits(f2, "try-error")) coef(f2)[3] else 0.25

    list(
      mu = c(0, 0),
      phi11 = max(min(phi11, 0.9), -0.9),
      phi12 = max(min(phi12, 0.9), -0.9),
      phi21 = max(min(phi21, 0.9), -0.9),
      phi22 = max(min(phi22, 0.9), -0.9),
      eta = rep(log(1), 2),
      rho_raw = rnorm(1, 0, 0.5)
    )
  }

  # ---- enumerate datasets to fit ---------------------------------------
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(paths)) {
    message("No simulation files found in ", data_dir)
    return(invisible())
  }

  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition & as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  if (!length(paths)) {
    message("Nothing to fit in requested range.")
    return(invisible())
  }

  message("Found ", length(paths), " datasets to process.")
  message("Overwrite fits: ", overwrite_fits)

  # ---- worker -----------------------------------------------------------
  fit_one_dataset <- function(p) {
    mm <- stringr::str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(mm[2])
    rid <- as.integer(mm[3])

    ds <- try(readRDS(p), silent = TRUE)
    if (inherits(ds, "try-error")) {
      return(sprintf("ERROR cond%03d rep%03d : readRDS failed", cid, rid))
    }

    y <- as.matrix(ds$data[, c("y1", "y2")])
    skew_dir <- as.integer(ds$skew_signs)
    sdat <- list(T = ds$T, y = y, skew_direction = skew_dir)

    msgs <- character()

    for (code in c("EI", "EL")) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))

      # Skip existing fits unless overwrite requested
      if (!overwrite_fits && file.exists(fit_path)) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
        next
      }

      # deterministic seeds
      seed_stan <- as.integer(SEED_BASE_STAN + cid * 1000L + rid + ifelse(code == "EL", 100000L, 0L))
      inits <- lapply(seq_len(chains), function(ch) {
        chain_seed <- as.integer(SEED_BASE_R + cid * 10000L + rid * 10L + ch + ifelse(code == "EL", 100000L, 0L))
        if (identical(code, "EI")) make_init_EI(y, chain_seed) else make_init_EL(y, chain_seed)
      })

      # sidecar for reproducibility
      try(
        saveRDS(
          list(
            meta = list(model = code, cid = cid, rid = rid, time = Sys.time()),
            init = inits,
            seed_stan = seed_stan
          ),
          file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
        ),
        silent = TRUE
      )

      # sampling
      warn_msgs <- character()
      fit <- withCallingHandlers(
        tryCatch(
          rstan::sampling(
            models[[code]],
            data = sdat,
            chains = chains,
            iter = iter,
            warmup = warmup,
            seed = seed_stan,
            init = inits,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            refresh = 0,
            cores = 1,
            pars = pars_keep,
            include = TRUE,
            save_warmup = FALSE
          ),
          error = function(e) e
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      if (inherits(fit, "stanfit") && length(fit@sim) > 0 && length(fit@sim$samples) > 0) {
        # divergences
        divs <- tryCatch(
          {
            sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
            if (is.list(sp) && length(sp) > 0) {
              sum(vapply(sp, function(x) as.integer(sum(x[, "divergent__"])), 0L))
            } else {
              NA_integer_
            }
          },
          error = function(e) NA_integer_
        )

        summ <- tryCatch(summary(fit)$summary, error = function(e) NULL)
        if (is.null(summ)) {
          out <- list(
            status = "summary_failed",
            error = "summary(fit) failed",
            model = code,
            condition_id = cid,
            rep_id = rid,
            n_div = divs,
            max_rhat = NA_real_,
            warnings = if (length(warn_msgs)) unique(warn_msgs) else NULL
          )
          saveRDS(out, fit_path)
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – summary_failed", code, cid, rid))
          next
        }

        out <- list(
          summary = summ,
          n_div = divs,
          max_rhat = suppressWarnings(max(summ[, "Rhat"], na.rm = TRUE)),
          model = code,
          condition_id = cid,
          rep_id = rid,
          status = "ok",
          warnings = if (length(warn_msgs)) unique(warn_msgs) else NULL
        )
        saveRDS(out, fit_path)
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : OK (div=%s)",
          code, cid, rid, ifelse(is.na(divs), "?", as.character(divs))
        ))
      } else if (inherits(fit, "error")) {
        out <- list(
          status = "error",
          error = conditionMessage(fit),
          model = code,
          condition_id = cid,
          rep_id = rid,
          n_div = NA_integer_,
          max_rhat = NA_real_,
          warnings = if (length(warn_msgs)) unique(warn_msgs) else NULL
        )
        saveRDS(out, fit_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – %s", code, cid, rid, conditionMessage(fit)))
      } else {
        out <- list(
          status = "unexpected",
          class = class(fit)[1],
          model = code,
          condition_id = cid,
          rep_id = rid
        )
        saveRDS(out, fit_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – unexpected", code, cid, rid))
      }
    }

    msgs
  }

  # ---- run --------------------------------------------------------------
  options(mc.cores = 1) # each fit uses one core internally

  use_fork <- (.Platform$OS.type != "windows") && cores_outer > 1
  if (use_fork) {
    message("Using mclapply (forking) with ", cores_outer, " workers.")
    logs <- parallel::mclapply(
      paths,
      fit_one_dataset,
      mc.cores = cores_outer,
      mc.preschedule = FALSE
    )
  } else {
    message("Using sequential fitting (Windows or cores_outer <= 1).")
    logs <- lapply(paths, fit_one_dataset)
  }

  # ---- append logs ------------------------------------------------------
  flat <- unlist(logs, use.names = FALSE)
  if (length(flat)) {
    con <- try(base::file(log_file, "a"), silent = TRUE)
    if (!inherits(con, "try-error")) {
      on.exit({ if (isOpen(con)) close(con) }, add = TRUE)
      try(writeLines(flat, con), silent = TRUE)
    } else {
      message("Logger fallback (", log_file, "):")
      for (r in flat) message(r)
    }
  }

  message(
    "\nFinished SEM fitting: ", length(paths), " datasets.\n",
    "Log → ", log_file, "\n",
    "Init sidecars → ", side_dir
  )
}
