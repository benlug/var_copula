###########################################################################
# fit_models.R – Study III: fit NG (Gaussian copula) & NC (Clayton copula)
###########################################################################

fit_var1_copula_models <- function(data_dir, fits_dir, stan_dir, results_dir,
                                   chains = 4, iter = 3000, warmup = 1500,
                                   adapt_delta = 0.98, max_treedepth = 13,
                                   cores_outer = 2,
                                   start_condition = 1, start_rep = 1) {
  suppressPackageStartupMessages({
    if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan required.")
    for (p in c("stringr", "digest", "doParallel", "foreach", "dplyr")) {
      if (!requireNamespace(p, quietly = TRUE)) stop(p, " required.")
    }
    library(rstan)
    library(stringr)
    library(digest)
    library(doParallel)
    library(foreach)
    library(dplyr)
  })

  SEED_BASE_R <- 7e5
  SEED_BASE_STAN <- 8e5
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  side_dir <- file.path(results_dir, "init_sidecar")
  dir.create(fits_dir, FALSE, TRUE)
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  # -- model cache by source hash ----------------------------------------
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
      message("Loading cached model: ", tag, " [", substr(h, 1, 8), "]")
      readRDS(rds)
    } else {
      message("Compiling model: ", tag, " [", substr(h, 1, 8), "]")
      m <- rstan::stan_model(path, model_name = paste0(tag, "_", substr(h, 1, 8)))
      saveRDS(m, rds)
      m
    }
  }

  models <- list(
    NG = model_cache(file.path(stan_dir, "model_NG_sl.stan"), "Normal_Gauss"),
    NC = model_cache(file.path(stan_dir, "model_NC_sl.stan"), "Normal_Clayton")
  )
  models <- Filter(Negate(is.null), models)
  if (!length(models)) stop("No models loaded. Check Stan files.")
  message("Successfully loaded models: ", paste(names(models), collapse = ", "))

  # -- init helpers -------------------------------------------------------
  make_init_common <- function() {
    list(
      mu    = rnorm(2, 0, 0.1),
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      sigma = runif(2, 0.8, 1.2)
    )
  }
  make_init <- function(code) {
    b <- make_init_common()
    if (code == "NG") {
      c(b, list(rho = runif(1, -0.5, 0.5)))
    } else if (code == "NC") {
      c(b, list(theta = runif(1, 0.2, 2.5)))
    } else {
      NULL
    }
  }

  # -- enumerate datasets & resume ---------------------------------------
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
    message("Nothing to fit (all requested cells done).")
    return(invisible())
  }

  # -- parallel outer-loop -----------------------------------------------
  if (is.null(cores_outer) || !is.numeric(cores_outer) || cores_outer < 1) cores_outer <- 1
  message("Starting parallel cluster with ", cores_outer, " cores.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  })

  options(mc.cores = 1) # avoid nested parallelism inside workers

  log_vec <- foreach::foreach(
    p = paths, .packages = c("rstan", "digest", "stringr", "dplyr"), .errorhandling = "pass"
  ) %dopar% {
    m <- stringr::str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(m[2])
    rid <- as.integer(m[3])

    ds <- try(readRDS(p), silent = TRUE)
    if (inherits(ds, "try-error")) {
      return(sprintf("ERROR cond%03d rep%03d : Failed to read data file.", cid, rid))
    }

    y_matrix <- as.matrix(ds$data[, c("y1", "y2")])
    sdat_base <- list(T = ds$T, y = y_matrix)

    msgs <- character()
    models_to_run <- intersect(c("NG", "NC"), names(models))
    for (code in models_to_run) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))

      if (file.exists(fit_path)) {
        chk <- try(readRDS(fit_path), silent = TRUE)
        if (inherits(chk, "stanfit") && length(chk@sim) > 0 && length(chk@sim$samples) > 0) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
          next
        } else {
          try(unlink(fit_path), silent = TRUE)
        }
      }

      set.seed(SEED_BASE_R + cid * 1000 + rid)
      inits <- replicate(chains, make_init(code), simplify = FALSE)
      if (!length(inits) || is.null(inits[[1]])) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – init generation failed", code, cid, rid))
        next
      }

      seeds_stan <- SEED_BASE_STAN + cid * 1000 + rid + 0:(chains - 1)
      try(saveRDS(
        list(
          meta = list(model = code, cid = cid, rid = rid, time = Sys.time()),
          seeds_stan = seeds_stan, init = inits
        ),
        file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
      ), silent = TRUE)

      warn_msgs <- character()
      fit <- withCallingHandlers(
        tryCatch(
          rstan::sampling(
            models[[code]],
            data = sdat_base,
            chains = chains, iter = iter, warmup = warmup, seed = seeds_stan[1],
            init = inits, control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            refresh = 0, cores = 1
          ),
          error = function(e) e
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      if (inherits(fit, "stanfit") && length(fit@sim) > 0 && length(fit@sim$samples) > 0) {
        saveRDS(fit, fit_path)
        divs <- tryCatch(
          {
            sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
            if (is.list(sp) && length(sp) > 0) sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L)) else NA_integer_
          },
          error = function(e) NA_integer_
        )
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : OK (%s div)", code, cid, rid, ifelse(is.na(divs), "NA", as.character(divs))))
        if (length(warn_msgs)) msgs <- c(msgs, sprintf("%s cond%03d rep%03d : WARN – %s", code, cid, rid, paste(unique(warn_msgs), collapse = " | ")))
      } else if (inherits(fit, "error")) {
        fail_path <- file.path(fits_dir, sprintf("fit_ERR_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – %s", code, cid, rid, conditionMessage(fit)))
      } else {
        fail_path <- file.path(fits_dir, sprintf("fit_MISC_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – unexpected return type", code, cid, rid))
      }
    }
    msgs
  }

  # -- write log ----------------------------------------------------------
  log_entries <- character()
  errors_occurred <- FALSE
  for (result in log_vec) {
    if (inherits(result, "error")) {
      log_entries <- c(log_entries, paste("FOREACH ERROR:", conditionMessage(result)))
      errors_occurred <- TRUE
    } else if (is.character(result)) {
      log_entries <- c(log_entries, result)
    }
  }
  if (length(log_entries)) {
    con <- base::file(log_file, "a")
    on.exit(
      {
        if (inherits(con, "connection")) {
          io <- try(isOpen(con), silent = TRUE)
          if (!inherits(io, "try-error") && io) close(con)
        }
      },
      add = TRUE
    )
    try(writeLines(log_entries, con), silent = TRUE)
  }

  message(
    "\nFinished processing ", length(paths), " data sets.",
    "\nLog      → ", log_file,
    "\nSide‑cars → ", side_dir
  )
  if (errors_occurred) message("!!! Errors occurred during parallel execution. Check the log file.")
}
