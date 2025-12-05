#!/usr/bin/env Rscript
###########################################################################
# fit_models.R  – robust Stan fitting (single‑level)
#   • unique seeds per condition × replication
#   • no nested parallelism
#   • distinct seeds for R and Stan
#   • efficient data input (matrix)
###########################################################################

fit_var1_copula_models <- function(data_dir,
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
                                   start_rep = 1) {
  suppressPackageStartupMessages({
    library(rstan)
    library(stringr)
    library(digest)
    library(doParallel)
    library(foreach)
    library(dplyr)
  })

  # FIX 3: Distinct seed bases for R (initialization) and Stan (sampling)
  SEED_BASE_R <- 1e6
  SEED_BASE_STAN <- 5e6

  side_dir <- file.path(results_dir, "init_sidecar")
  dir.create(fits_dir, FALSE, TRUE)
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  ## -- compile or load cached models ------------------------------------
  model_cache <- function(path, tag) {
    rds <- file.path(stan_dir, paste0(tag, ".rds"))
    if (file.exists(rds)) {
      readRDS(rds)
    } else {
      # NOTE: Ensure Stan models are updated to accept data 'matrix[T, 2] y;'
      # instead of 'vector[2] y[T];' to utilize Fix 5.
      m <- stan_model(path, model_name = tag)
      saveRDS(m, rds)
      m
    }
  }
  models <- list(
    SG = model_cache(file.path(stan_dir, "model_SNG_sl.stan"), "SkewN_Gauss"),
    NG = model_cache(file.path(stan_dir, "model_NG_sl.stan"), "Normal_Gauss")
  )

  ## -- helper: minimal init list per model ------------------------------
  # Improvement: Initialize closer to expected values to aid fitting.
  # In fit_models.R, replace the make_init function:

  ## -- helper: minimal init list per model ------------------------------
  # CRITICAL FIX: Initialize parameters matching the Stan model definitions.
  make_init <- function(code) {
    # Base parameters common to both models (mu, phi, rho)
    base <- list(
      mu = rnorm(2, 0, 0.1),
      # Initialize phi within a reasonably stationary region
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      rho = runif(1, -0.5, 0.5)
    )

    # NOTE: 'sigma' must be removed from 'base' as it is only used in the NG model.

    if (code == "SG") {
      # SG Model (CP): Initialize omega and delta.
      c(base, list(
        omega = runif(2, 0.8, 1.5), # Omega > 0
        # delta is bounded (-1, 1). Initialize away from boundaries.
        delta = runif(2, -0.8, 0.8)
      ))
    } else { # code == "NG"
      # NG Model: Initialize sigma.
      c(base, list(
        # Initialize sigma near 1 as data is standardized
        sigma = runif(2, 0.8, 1.2)
      ))
    }
  }

  ## -- list data sets & honour resume -----------------------------------
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  meta <- str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition &
      as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  if (!length(paths)) {
    message("Nothing to fit (all requested cells done).")
    return(invisible())
  }

  ## -- foreach cluster ---------------------------------------------------
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    parallel::stopCluster(cl)
    doParallel::registerDoSEQ()
  })

  options(mc.cores = 1) # Stan uses ONE core

  log_vec <- foreach(
    p              = paths,
    .packages      = c("rstan", "digest", "stringr"),
    .errorhandling = "pass"
  ) %dopar% {
    m <- str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(m[2])
    rid <- as.integer(m[3])

    ds <- readRDS(p)

    # FIX 5: Prepare data as a T x 2 matrix instead of a list of vectors
    sdat <- list(
      T = ds$T,
      y = as.matrix(ds$data[, c("y1", "y2")])
    )

    msgs <- character()

    for (code in names(models)) {
      fit_path <- file.path(
        fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid)
      )
      if (file.exists(fit_path)) {
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : skip (exists)",
          code, cid, rid
        ))
        next
      }

      ## --- robust init & unique seeds ----------------------------------
      # FIX 3: Use SEED_BASE_R for R operations (generating inits)
      seed_r <- SEED_BASE_R + cid * 1000 + rid
      set.seed(seed_r)
      inits <- replicate(chains, make_init(code), simplify = FALSE)

      # FIX 3: Use SEED_BASE_STAN for Stan sampling
      seeds_stan <- SEED_BASE_STAN + cid * 1000 + rid + 0:(chains - 1)

      saveRDS(
        list(
          meta = list(
            model = code, condition_id = cid,
            rep_id = rid, time = Sys.time()
          ),
          seeds_r = seed_r,
          seeds_stan = seeds_stan,
          init = inits
        ),
        file.path(
          side_dir,
          sprintf(
            "INIT_%s_cond%03d_rep%03d.rds",
            code, cid, rid
          )
        )
      )

      ## --- sampling ----------------------------------------------------
      fit <- tryCatch(
        sampling(models[[code]],
          data = sdat,
          chains = chains, iter = iter, warmup = warmup,
          # Stan uses the first seed; others are derived internally if cores > 1
          seed = seeds_stan[1],
          init = inits,
          control = list(
            adapt_delta = adapt_delta,
            max_treedepth = max_treedepth
          ),
          refresh = 0,
          cores = 1
        ), # enforce single‑core
        error = function(e) e
      )

      saveRDS(fit, fit_path)

      if (inherits(fit, "stanfit")) {
        divs <- sum(vapply(
          get_sampler_params(fit, FALSE),
          \(x) sum(x[, "divergent__"]), 0
        ))
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : OK (%d div)",
          code, cid, rid, divs
        ))
      } else {
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : FAIL – %s",
          code, cid, rid, conditionMessage(fit)
        ))
      }
    }
    msgs
  }

  ## -- write log safely --------------------------------------------------
  # Appending to the log file for safety across resumes/crashes
  if (length(log_vec) > 0) {
    con <- file(log_file, "a")
    writeLines(unlist(log_vec), con)
    close(con)
  }

  message(
    "\nFinished ", length(paths), " data sets.",
    "\nLog      → ", log_file,
    "\nSide‑cars → ", side_dir
  )
}
