#!/usr/bin/env Rscript
###########################################################################
# fit_models.R  - Robust Stan fitting (single-level)
#   - Compiles models INSIDE workers to avoid serialization errors
#   - Saves only essential fit summaries (not full stanfit objects)
#   - FIXED: Stan cache staleness check using source file hash
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
    library(stringr)
    library(doParallel)
    library(foreach)
    library(dplyr)
  })

  # Seed bases for reproducibility
  SEED_BASE_R <- 1e6
  SEED_BASE_STAN <- 5e6

  dir.create(fits_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  ## -- List data sets & honor resume --------------------------------------
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$",
                      full.names = TRUE)
  meta <- str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition &
       as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  
  if (!length(paths)) {
    message("Nothing to fit (all requested cells done).")
    return(invisible())
  }
  
  message("Fitting ", length(paths), " data sets using ", cores_outer, " cores...")

  ## -- Foreach cluster ----------------------------------------------------
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    parallel::stopCluster(cl)
    doParallel::registerDoSEQ()
  })

  # Export necessary variables to workers (but NOT compiled models)
  parallel::clusterExport(cl, c("stan_dir", "fits_dir", "chains", "iter", 
                                 "warmup", "adapt_delta", "max_treedepth",
                                 "SEED_BASE_R", "SEED_BASE_STAN"),
                          envir = environment())

  log_vec <- foreach(
    p = paths,
    .packages = c("rstan", "stringr", "dplyr", "tibble"),
    .errorhandling = "pass"
  ) %dopar% {
    
    # ── Helper: compute file hash for cache validation ──
    file_hash <- function(path) {
      if (!file.exists(path)) return(NULL)
      # Use md5sum if digest not available
      if (requireNamespace("digest", quietly = TRUE)) {
        digest::digest(file = path, algo = "md5")
      } else {
        # Fallback: use file modification time + size
        info <- file.info(path)
        paste0(info$mtime, "_", info$size)
      }
    }
    
    # ── Helper: compile or load cached model with staleness check ──
    model_cache <- function(stan_path, tag) {
      rds_path <- file.path(stan_dir, paste0(tag, ".rds"))
      hash_path <- file.path(stan_dir, paste0(tag, ".hash"))
      
      current_hash <- file_hash(stan_path)
      
      # Check if cache exists and is valid
      cache_valid <- FALSE
      if (file.exists(rds_path) && file.exists(hash_path)) {
        cached_hash <- readLines(hash_path, n = 1, warn = FALSE)
        cache_valid <- identical(cached_hash, current_hash)
      }
      
      if (cache_valid) {
        return(readRDS(rds_path))
      }
      
      # Compile and cache
      m <- rstan::stan_model(stan_path, model_name = tag)
      saveRDS(m, rds_path)
      writeLines(current_hash, hash_path)
      m
    }
    
    # ── Helper: initialization list per model ──
    make_init <- function(code) {
      base <- list(
        mu = rnorm(2, 0, 0.1),
        phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
        phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
        rho = runif(1, -0.5, 0.5)
      )
      if (code == "SG") {
        c(base, list(omega = runif(2, 0.8, 1.5), delta = runif(2, -0.8, 0.8)))
      } else {
        c(base, list(sigma = runif(2, 0.8, 1.2)))
      }
    }
    
    # ── Helper: extract essential information from stanfit ──
    extract_fit_summary <- function(fit, code, cid, rid) {
      if (!inherits(fit, "stanfit")) {
        return(list(
          meta = list(model = code, condition_id = cid, rep_id = rid,
                      status = "error", 
                      error_msg = if(inherits(fit, "error")) conditionMessage(fit) else "unknown"),
          summary = NULL,
          diagnostics = NULL
        ))
      }
      
      # Check if sampling actually occurred
      if (length(fit@sim) == 0 || is.null(fit@sim$iter) || 
          fit@sim$iter == 0 || length(fit@sim$samples) == 0) {
        return(list(
          meta = list(model = code, condition_id = cid, rep_id = rid, status = "empty"),
          summary = NULL,
          diagnostics = NULL
        ))
      }
      
      # Extract parameter summary
      s <- tryCatch(rstan::summary(fit)$summary, error = function(e) NULL)
      
      if (is.null(s)) {
        return(list(
          meta = list(model = code, condition_id = cid, rep_id = rid, status = "summary_error"),
          summary = NULL,
          diagnostics = NULL
        ))
      }
      
      # Extract diagnostics
      sp <- tryCatch(rstan::get_sampler_params(fit, FALSE), error = function(e) NULL)
      n_div <- if (!is.null(sp)) {
        sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L))
      } else NA_integer_
      
      max_rhat <- suppressWarnings(max(s[, "Rhat"], na.rm = TRUE))
      
      # Extract ESS (new diagnostic)
      min_ess_bulk <- suppressWarnings(min(s[, "n_eff"], na.rm = TRUE))
      
      # Convert summary to data frame
      summary_df <- as.data.frame(s) |>
        tibble::rownames_to_column("param") |>
        dplyr::select(param, mean, sd, `2.5%`, `25%`, `50%`, `75%`, `97.5%`, n_eff, Rhat)
      
      list(
        meta = list(
          model = code, condition_id = cid, rep_id = rid,
          status = "ok", time = Sys.time()
        ),
        summary = summary_df,
        diagnostics = list(
          n_div = n_div, 
          max_rhat = max_rhat,
          min_ess = min_ess_bulk,
          n_chains = fit@sim$chains, 
          iter = fit@sim$iter, 
          warmup = fit@sim$warmup
        )
      )
    }
    
    # ── Main worker logic ──
    m <- str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(m[2])
    rid <- as.integer(m[3])

    ds <- readRDS(p)
    sdat <- list(T = ds$T, y = as.matrix(ds$data[, c("y1", "y2")]))

    msgs <- character()
    
    # Load models inside worker (avoids serialization issues)
    models <- list(
      SG = model_cache(file.path(stan_dir, "model_SNG_sl.stan"), "SkewN_Gauss"),
      NG = model_cache(file.path(stan_dir, "model_NG_sl.stan"), "Normal_Gauss")
    )

    for (code in names(models)) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))
      
      if (file.exists(fit_path)) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
        next
      }

      # Generate inits with reproducible seed
      seed_r <- SEED_BASE_R + cid * 1000 + rid
      set.seed(seed_r)
      inits <- replicate(chains, make_init(code), simplify = FALSE)

      # Stan sampling seed
      seed_stan <- SEED_BASE_STAN + cid * 1000 + rid

      # Sampling
      fit <- tryCatch(
        rstan::sampling(models[[code]],
                        data = sdat,
                        chains = chains, iter = iter, warmup = warmup,
                        seed = seed_stan,
                        init = inits,
                        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                        refresh = 0,
                        cores = 1),
        error = function(e) e
      )

      # Extract and save only essential information
      fit_summary <- extract_fit_summary(fit, code, cid, rid)
      saveRDS(fit_summary, fit_path)

      if (fit_summary$meta$status == "ok") {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : OK (%d div, ESS=%.0f)", 
                                code, cid, rid, 
                                fit_summary$diagnostics$n_div,
                                fit_summary$diagnostics$min_ess))
      } else {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : %s", 
                                code, cid, rid, fit_summary$meta$status))
      }
    }
    msgs
  }

  ## -- Write log ----------------------------------------------------------
  if (length(log_vec) > 0) {
    con <- file(log_file, "a")
    writeLines(unlist(log_vec), con)
    close(con)
  }

  message("\nFinished ", length(paths), " data sets.\nLog -> ", log_file)
}
