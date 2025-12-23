#!/usr/bin/env Rscript
###########################################################################
# fit_models.R  – EXACT COPY OF TVP STRUCTURE
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
    if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan package required.")
    if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr package required.")
    if (!requireNamespace("digest", quietly = TRUE)) stop("digest package required.")
    if (!requireNamespace("doParallel", quietly = TRUE)) stop("doParallel package required.")
    if (!requireNamespace("foreach", quietly = TRUE)) stop("foreach package required.")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr package required.")

    library(rstan)
    library(stringr)
    library(digest)
    library(doParallel)
    library(foreach)
    library(dplyr)
  })

  SEED_BASE_R <- 1e6
  SEED_BASE_STAN <- 5e6
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  side_dir <- file.path(results_dir, "init_sidecar")
  dir.create(fits_dir, FALSE, TRUE)
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  ## -- model cache (hash by source) ------------------------------------------
  model_cache <- function(path, tag) {
    if (!file.exists(path)) {
      message("Model file not found: ", path)
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

  ## -- load models -----------------------------------------------------------
  models <- list(
    SG = model_cache(file.path(stan_dir, "model_SNG_sl.stan"), "SkewN_Gauss"),
    NG = model_cache(file.path(stan_dir, "model_NG_sl.stan"), "Normal_Gauss")
  )
  models <- Filter(Negate(is.null), models)
  
  if (!length(models)) stop("No models loaded successfully. Check Stan files.")
  message("Successfully loaded models: ", paste(names(models), collapse = ", "))

  ## -- init functions --------------------------------------------------------
  make_init_SG <- function() {
    list(
      mu = rnorm(2, 0, 0.1),
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      omega = runif(2, 0.8, 1.5), 
      delta = runif(2, -0.8, 0.8),
      rho = runif(1, -0.5, 0.5)
    )
  }
  
  make_init_NG <- function() {
    list(
      mu = rnorm(2, 0, 0.1),
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      sigma = runif(2, 0.8, 1.2),
      rho = runif(1, -0.5, 0.5)
    )
  }

  ## -- enumerate datasets ----------------------------------------------------
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

  ## -- parallel setup --------------------------------------------------------
  if (is.null(cores_outer) || !is.numeric(cores_outer) || cores_outer < 1) cores_outer <- 1
  message("Starting parallel cluster with ", cores_outer, " cores.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    if (requireNamespace("foreach", quietly = TRUE)) try(foreach::registerDoSEQ(), silent = TRUE)
  })

  options(mc.cores = 1)

  ## -- main fitting loop -----------------------------------------------------
  log_vec <- foreach::foreach(
    p = paths,
    .packages = c("rstan", "digest", "stringr", "dplyr"),
    .errorhandling = "pass"
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

    for (code in names(models)) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))

      # Skip if exists and valid (compact summary format)
      if (file.exists(fit_path)) {
        chk <- try(readRDS(fit_path), silent = TRUE)
        if (is.list(chk) && !is.null(chk$status) && chk$status == "ok" && !is.null(chk$summary)) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
          next
        } else {
          try(unlink(fit_path), silent = TRUE)
        }
      }

      # Prepare data and inits
      sdat <- sdat_base
      
      if (code == "SG") {
        inits <- replicate(chains, make_init_SG(), simplify = FALSE)
      } else if (code == "NG") {
        inits <- replicate(chains, make_init_NG(), simplify = FALSE)
      }

      if (!length(inits) || is.null(inits[[1]])) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – init generation failed", code, cid, rid))
        next
      }

      seeds_stan <- SEED_BASE_STAN + cid * 1000 + rid + 0:(chains - 1)

      # Save init sidecar
      try(saveRDS(
        list(
          meta = list(model = code, condition_id = cid, rep_id = rid, time = Sys.time()),
          seeds_stan = seeds_stan, init = inits
        ),
        file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
      ), silent = TRUE)

      ## --- sampling ---------------------------------------------------------
      warn_msgs <- character()
      fit <- withCallingHandlers(
        tryCatch(
          rstan::sampling(
            models[[code]],
            data = sdat,
            chains = chains, iter = iter, warmup = warmup,
            seed = seeds_stan[1],
            init = inits,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            refresh = 0,
            cores = 1
          ),
          error = function(e) e
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      if (inherits(fit, "stanfit") && length(fit@sim) > 0 && length(fit@sim$samples) > 0) {
        
        # Extract only what we need for analysis (much smaller than full fit)
        divs <- tryCatch({
          sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
          if (is.list(sp) && length(sp) > 0) sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L)) else NA_integer_
        }, error = function(e) NA_integer_)
        
        max_rhat_val <- tryCatch({
          s <- summary(fit)$summary
          suppressWarnings(max(s[, "Rhat"], na.rm = TRUE))
        }, error = function(e) NA_real_)
        
        # Save compact summary object instead of full fit
        fit_summary <- tryCatch({
          list(
            summary = summary(fit)$summary,
            n_div = divs,
            max_rhat = max_rhat_val,
            sampler_params = list(
              n_divergent = divs,
              max_treedepth_hits = tryCatch({
                sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
                sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= max_treedepth), 0L))
              }, error = function(e) NA_integer_)
            ),
            model = code,
            condition_id = cid,
            rep_id = rid,
            status = "ok"
          )
        }, error = function(e) {
          list(status = "summary_failed", error = conditionMessage(e))
        })
        
        saveRDS(fit_summary, fit_path)

        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : OK (%s div)", code, cid, rid,
                                ifelse(is.na(divs), "NA", as.character(divs))))
        if (length(warn_msgs)) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : WARN – %s", code, cid, rid,
                                  paste(unique(warn_msgs), collapse = " | ")))
        }
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

  ## -- write log -------------------------------------------------------------
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
    on.exit({
      if (inherits(con, "connection")) {
        io <- try(isOpen(con), silent = TRUE)
        if (!inherits(io, "try-error") && io) close(con)
      }
    }, add = TRUE)
    try(writeLines(log_entries, con), silent = TRUE)
  }

  message("\nFinished processing ", length(paths), " data sets.",
          "\nLog      → ", log_file,
          "\nSide-cars → ", side_dir)
  if (errors_occurred) {
    message("!!! Errors occurred during parallel execution. Check the log file.")
  }
}
