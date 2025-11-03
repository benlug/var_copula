#!/usr/bin/env Rscript
###########################################################################
# fit_models_SEM.R — Stan fitting for SEM Study A/B
# Codes:
#   EI = Exponential margins at Indicator layer (Study A)
#   EL = Exponential margins at Latent/State layer (Study B, ε ≡ 0)
# Mirrors Study 2 compilation/caching & logging. :contentReference[oaicite:4]{index=4}
###########################################################################

fit_sem_models <- function(data_dir, fits_dir, stan_dir, results_dir,
                           chains = 4, iter = 3000, warmup = 1500,
                           adapt_delta = 0.995, max_treedepth = 15,
                           cores_outer = 2,
                           start_condition = 1, start_rep = 1) {
  suppressPackageStartupMessages({
    pkgs <- c("rstan", "stringr", "digest", "doParallel", "foreach", "dplyr")
    for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) stop(p, " required.")
    library(rstan)
    library(stringr)
    library(digest)
    library(doParallel)
    library(foreach)
    library(dplyr)
  })

  SEED_BASE_R <- 7e6
  SEED_BASE_STAN <- 9e6
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  side_dir <- file.path(results_dir, "init_sidecar_sem")
  dir.create(fits_dir, FALSE, TRUE)
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit_sem.log")

  # ---- compile/load cached models by hash -------------------------------
  model_cache <- function(path, tag) {
    if (!file.exists(path)) {
      return(NULL)
    }
    code <- readChar(path, file.info(path)$size)
    h <- digest::digest(code, algo = "xxhash64")
    rds <- file.path(stan_dir, sprintf("%s_%s.rds", tag, substr(h, 1, 8)))
    old <- list.files(stan_dir, paste0("^", tag, "_[0-9a-f]{8}\\.rds$"), full.names = TRUE)
    old <- setdiff(old, rds)
    if (length(old)) try(unlink(old), silent = TRUE)
    if (file.exists(rds)) {
      message("Loading cached: ", tag, " [", substr(h, 1, 8), "]")
      readRDS(rds)
    } else {
      message("Compiling: ", tag, " [", substr(h, 1, 8), "]")
      m <- rstan::stan_model(path, model_name = paste0(tag, "_", substr(h, 1, 8)))
      saveRDS(m, rds)
      m
    }
  }

  models <- list(
    EI = model_cache(file.path(stan_dir, "sem_A_indicator_EG.stan"), "SEM_Indicator_EG"),
    EL = model_cache(file.path(stan_dir, "sem_B_latent_EG.stan"), "SEM_Latent_EG")
  )
  models <- Filter(Negate(is.null), models)
  if (!length(models)) stop("No SEM models loaded. Check stan files.")

  # ---- init helpers -----------------------------------------------------
  make_init_EI <- function(y, chain_id) {
    set.seed(SEED_BASE_R + chain_id)
    list(
      mu = rnorm(2, 0, 0.1),
      phi11 = 0.55, phi12 = 0.10, phi21 = 0.10, phi22 = 0.25,
      eta = rep(log(1.0), 2),
      rho = runif(1, -0.2, 0.2),
      # provide a rough latent trajectory to stabilize the Gaussian state equation
      state = as.matrix(y)
    )
  }
  make_init_EL <- function(y, chain_id) {
    set.seed(SEED_BASE_R + chain_id)
    # rough OLS to seed phi
    ylag <- rbind(c(0, 0), y[-nrow(y), , drop = FALSE])
    f1 <- try(lm(y[, 1] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
    f2 <- try(lm(y[, 2] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
    phi11 <- if (!inherits(f1, "try-error")) coef(f1)[2] %||% 0.55 else 0.55
    phi12 <- if (!inherits(f1, "try-error")) coef(f1)[3] %||% 0.10 else 0.10
    phi21 <- if (!inherits(f2, "try-error")) coef(f2)[2] %||% 0.10 else 0.10
    phi22 <- if (!inherits(f2, "try-error")) coef(f2)[3] %||% 0.25 else 0.25
    list(
      mu = c(0, 0),
      phi11 = phi11, phi12 = phi12, phi21 = phi21, phi22 = phi22,
      eta = rep(log(1.0), 2),
      rho = runif(1, -0.2, 0.2)
    )
  }

  # ---- enumerate datasets & resume -------------------------------------
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  if (!length(paths)) {
    message("No data files found in ", data_dir)
    return(invisible())
  }

  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition & as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  if (!length(paths)) {
    message("Nothing to fit in requested range.")
    return(invisible())
  }

  # ---- outer parallel ---------------------------------------------------
  cores_outer <- if (is.null(cores_outer) || cores_outer < 1) 1 else as.integer(cores_outer)
  message("Starting parallel cluster with ", cores_outer, " cores.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    if (requireNamespace("foreach", quietly = TRUE)) try(foreach::registerDoSEQ(), silent = TRUE)
  })

  options(mc.cores = 1) # avoid nested inside workers

  log_vec <- foreach::foreach(
    p = paths, .packages = c("rstan", "digest", "stringr", "dplyr"), .errorhandling = "pass"
  ) %dopar% {
    mm <- stringr::str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(mm[2])
    rid <- as.integer(mm[3])
    ds <- try(readRDS(p), silent = TRUE)
    if (inherits(ds, "try-error")) {
      return(sprintf("ERROR cond%03d rep%03d: readRDS failed", cid, rid))
    }

    y <- as.matrix(ds$data[, c("y1", "y2")])
    sd <- list(T = ds$T, y = y, skew_direction = as.integer(ds$skew_signs))
    msgs <- character()

    for (code in c("EI", "EL")) {
      fpath <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))
      if (file.exists(fpath)) {
        chk <- try(readRDS(fpath), silent = TRUE)
        if (inherits(chk, "stanfit") && length(chk@sim) > 0 && length(chk@sim$samples) > 0) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
          next
        } else {
          try(unlink(fpath), silent = TRUE)
        }
      }

      inits <- switch(code,
        EI = lapply(seq_len(chains), function(k) make_init_EI(y, cid * 1000 + rid * 10 + k)),
        EL = lapply(seq_len(chains), function(k) make_init_EL(y, cid * 1000 + rid * 10 + k))
      )

      seeds <- SEED_BASE_STAN + cid * 1000 + rid + 0:(chains - 1)

      try(
        saveRDS(
          list(
            meta = list(model = code, cid = cid, rid = rid, time = Sys.time()),
            init = inits, seeds = seeds
          ),
          file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
        ),
        silent = TRUE
      )

      fit <- withCallingHandlers(
        tryCatch(
          rstan::sampling(models[[code]],
            data = sd,
            chains = chains, iter = iter, warmup = warmup, seed = seeds[1],
            init = inits,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            refresh = 0, cores = 1
          ),
          error = function(e) e
        ),
        warning = function(w) invokeRestart("muffleWarning")
      )

      if (inherits(fit, "stanfit") && length(fit@sim) > 0) {
        saveRDS(fit, fpath)
        divs <- tryCatch(
          {
            sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
            if (is.list(sp) && length(sp) > 0) sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L)) else NA_integer_
          },
          error = function(e) NA_integer_
        )
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : OK (%s div)",
          code, cid, rid, ifelse(is.na(divs), "NA", as.character(divs))
        ))
      } else {
        fail <- file.path(fits_dir, sprintf("fit_ERR_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL", code, cid, rid))
      }
    }
    msgs
  }

  if (length(log_vec)) {
    con <- base::file(log_file, "a")
    on.exit(
      {
        if (isOpen(con)) close(con)
      },
      add = TRUE
    )
    for (r in log_vec) if (is.character(r)) try(writeLines(r, con), silent = TRUE)
  }
  message("\nFinished SEM: ", length(paths), " datasets.\nLog → ", log_file, "\nSide‑cars → ", side_dir)
}
