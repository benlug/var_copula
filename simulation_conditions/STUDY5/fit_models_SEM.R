# fit_models_SEM.R — Fit Study 5 models with feasibility-lift Stan
# CORRECTED VERSION — Bug fixes applied

fit_sem_models <- function(data_dir, fits_dir, stan_dir, results_dir,
                           chains = 4, iter = 3000, warmup = 1500,
                           adapt_delta = 0.90, max_treedepth = 10,
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
  rstan::rstan_options(auto_write = TRUE)

  SEED_R <- 7e6
  SEED_STAN <- 9e6

  ## dirs & log
  side_dir <- file.path(results_dir, "init_sidecar_sem")
  dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(side_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(fits_dir, "stan_fit_sem.log")
  if (!file.exists(log_file)) try(file.create(log_file), silent = TRUE)

  ## compile / cache
  model_cache <- function(path, tag) {
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

  ## init helpers (feasible)
  make_init_EI <- function(y, skew_dir, chain_id) {
    set.seed(SEED_R + chain_id)
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
      eta = rep(log(0.2), 2), # modest slack
      rho_raw = rnorm(1, 0, 0.5),
      z_raw = matrix(0, nrow(y), 2)
    )
  }

  # ============================================================================
  # FIX: Corrected make_init_EL function
  # Bug 1: Was using 'rho' instead of 'rho_raw' (Stan parameter name mismatch)
  # Bug 2: Was missing 'eta' parameter entirely
  # ============================================================================
  make_init_EL <- function(y, skew_dir, chain_id) {
    set.seed(SEED_R + chain_id)
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
      eta = rep(log(1), 2),   # give more slack to keep sigma_exp positive
      rho_raw = rnorm(1, 0, 0.5) # FIX: Renamed from 'rho' to 'rho_raw'
    )
  }


  ## enumerate data
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(paths)) {
    message("No data files in ", data_dir)
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

  ## outer parallel
  message("Starting parallel cluster with ", cores_outer, " cores.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  })

  options(mc.cores = 1)

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
    skew_dir <- as.integer(ds$skew_signs)
    sd <- list(T = ds$T, y = y, skew_direction = skew_dir)
    msgs <- character()

    for (code in c("EI", "EL")) {
      fpath <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))
      errpth <- file.path(fits_dir, sprintf("fit_ERR_%s_cond%03d_rep%03d.rds", code, cid, rid))

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
        EI = lapply(seq_len(chains), function(k) make_init_EI(y, skew_dir, cid * 1000 + rid * 10 + k)),
        EL = lapply(seq_len(chains), function(k) make_init_EL(y, skew_dir, cid * 1000 + rid * 10 + k))
      )
      seeds <- SEED_STAN + cid * 1000 + rid + 0:(chains - 1)

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
            chains = chains, iter = iter, warmup = warmup,
            seed = seeds[1], init = inits,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            refresh = 0, cores = 1
          ),
          error = function(e) e
        ),
        warning = function(w) invokeRestart("muffleWarning")
      )

      if (inherits(fit, "stanfit") && length(fit@sim) > 0) {
        saveRDS(fit, fpath)
        if (file.exists(errpth)) try(unlink(errpth), silent = TRUE)
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
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : OK (%s div)",
          code, cid, rid, ifelse(is.na(divs), "NA", as.character(divs))
        ))
      } else {
        saveRDS(fit, errpth)
        msg <- if (inherits(fit, "error")) conditionMessage(fit) else paste("unexpected type:", class(fit)[1])
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – %s", code, cid, rid, msg))
      }
    }
    msgs
  }

  if (length(log_vec)) {
    con <- try(base::file(log_file, "a"), silent = TRUE)
    if (!inherits(con, "try-error")) {
      on.exit(
        {
          if (isOpen(con)) close(con)
        },
        add = TRUE
      )
      for (r in log_vec) if (is.character(r)) try(writeLines(r, con), silent = TRUE)
    } else {
      message("Logger fallback (", log_file, "):")
      for (r in log_vec) if (is.character(r)) message(r)
    }
  }

  message(
    "\nFinished SEM: ", length(paths), " datasets.\nLog → ", log_file,
    "\nSide-cars → ", file.path(results_dir, "init_sidecar_sem")
  )
}
