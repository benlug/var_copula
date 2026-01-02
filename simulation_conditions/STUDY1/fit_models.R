#!/usr/bin/env Rscript
###########################################################################
# fit_models.R  – Study 1: Using mclapply (forking) to avoid serialization
#
# mclapply uses forking which shares the parent's memory space, avoiding
# the serialization issues that corrupt Stan models with PSOCK clusters.
#
# NOTE: This works on Linux/Mac but NOT on Windows (Windows doesn't support fork)
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
    library(parallel)
    library(dplyr)
  })

  SEED_BASE_R <- 1e6
  SEED_BASE_STAN <- 5e6

  side_dir <- file.path(results_dir, "init_sidecar")
  dir.create(fits_dir, FALSE, TRUE)
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  ## -- Compile models in main process ----------------------------------------
  message("Compiling Stan models...")

  sg_path <- file.path(stan_dir, "model_SNG_sl.stan")
  ng_path <- file.path(stan_dir, "model_NG_sl.stan")

  if (!file.exists(sg_path)) stop("SG model not found: ", sg_path)
  if (!file.exists(ng_path)) stop("NG model not found: ", ng_path)

  model_SG <- rstan::stan_model(sg_path, model_name = "SkewN_Gauss")
  message("  Compiled: SG")

  model_NG <- rstan::stan_model(ng_path, model_name = "Normal_Gauss")
  message("  Compiled: NG")

  # Store in environment so forked processes can access them
  models <- list(SG = model_SG, NG = model_NG)

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
    message("Nothing to fit.")
    return(invisible())
  }

  message("Found ", length(paths), " datasets to process.")
  message("Using ", cores_outer, " cores with mclapply (forking).")

  ## -- Define worker function ------------------------------------------------
  fit_one_dataset <- function(p) {
    m <- stringr::str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(m[2])
    rid <- as.integer(m[3])

    ds <- try(readRDS(p), silent = TRUE)
    if (inherits(ds, "try-error")) {
      return(sprintf("ERROR cond%03d rep%03d : Failed to read data file.", cid, rid))
    }

    y_matrix <- as.matrix(ds$data[, c("y1", "y2")])
    sdat <- list(T = ds$T, y = y_matrix)

    msgs <- character()

    for (code in names(models)) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))

      # Skip if exists and valid
      if (file.exists(fit_path)) {
        chk <- try(readRDS(fit_path), silent = TRUE)
        if (is.list(chk) && !is.null(chk$status) && chk$status == "ok" && !is.null(chk$summary)) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip", code, cid, rid))
          next
        } else {
          try(unlink(fit_path), silent = TRUE)
        }
      }

      # Generate inits
      inits <- lapply(seq_len(chains), function(ch) {
        set.seed(SEED_BASE_R + cid * 10000L + rid * 10L + ch)
        if (identical(code, "SG")) make_init_SG() else make_init_NG()
      })

      seed_stan <- SEED_BASE_STAN + cid * 1000 + rid

      # Save init sidecar
      try(saveRDS(
        list(meta = list(model = code, cid = cid, rid = rid), init = inits),
        file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
      ), silent = TRUE)

      ## --- sampling (models available via forking) --------------------------
      warn_msgs <- character()
      fit <- withCallingHandlers(
        tryCatch(
          rstan::sampling(
            models[[code]], # Available because mclapply forks (shares memory)
            data = sdat,
            chains = chains, iter = iter, warmup = warmup,
            thin = 2,
            seed = seed_stan,
            init = inits,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            refresh = 0,
            cores = 1,
            save_warmup = FALSE
          ),
          error = function(e) e
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      ## --- save results -----------------------------------------------------
      if (inherits(fit, "stanfit") && length(fit@sim) > 0) {
        # Extract divergences - try multiple methods
        divs <- NA_integer_
        treedepth_hits <- NA_integer_

        # Method 1: Standard sampler params
        tryCatch(
          {
            sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
            if (is.list(sp) && length(sp) > 0) {
              # Check if first chain has the expected structure
              if (is.matrix(sp[[1]]) && "divergent__" %in% colnames(sp[[1]])) {
                divs <- as.integer(sum(sapply(sp, function(x) sum(x[, "divergent__"]))))
              }
              if (is.matrix(sp[[1]]) && "treedepth__" %in% colnames(sp[[1]])) {
                treedepth_hits <- as.integer(sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= max_treedepth))))
              }
            }
          },
          error = function(e) NULL
        )

        # Method 2: If still NA, try parsing from warnings (Stan warns about divergences)
        if (is.na(divs) && length(warn_msgs) > 0) {
          div_warn <- grep("divergent transition", warn_msgs, value = TRUE, ignore.case = TRUE)
          if (length(div_warn) > 0) {
            # Extract number from warning like "X of Y iterations..."
            nums <- as.integer(gsub(".*?([0-9]+).*", "\\1", div_warn[1]))
            if (!is.na(nums[1])) divs <- nums[1]
          }
        }

        # If no divergence info found, assume 0 (Stan warns if there are any)
        if (is.na(divs) && !any(grepl("divergent", warn_msgs, ignore.case = TRUE))) {
          divs <- 0L
        }

        fit_summary <- tryCatch(
          {
            summ <- summary(fit)$summary
            list(
              summary = summ,
              n_div = divs,
              max_rhat = suppressWarnings(max(summ[, "Rhat"], na.rm = TRUE)),
              sampler_params = list(
                n_divergent = divs,
                max_treedepth_hits = treedepth_hits
              ),
              warnings = if (length(warn_msgs)) unique(warn_msgs) else NULL,
              model = code,
              condition_id = cid,
              rep_id = rid,
              status = "ok"
            )
          },
          error = function(e) {
            list(
              status = "summary_failed", error = conditionMessage(e),
              model = code, condition_id = cid, rep_id = rid
            )
          }
        )

        saveRDS(fit_summary, fit_path)
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : OK (div=%s, tree=%s)", code, cid, rid,
          ifelse(is.na(divs), "?", as.character(divs)),
          ifelse(is.na(treedepth_hits), "?", as.character(treedepth_hits))
        ))
      } else if (inherits(fit, "error")) {
        saveRDS(list(
          status = "error", error = conditionMessage(fit),
          model = code, condition_id = cid, rep_id = rid
        ), fit_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – %s", code, cid, rid, conditionMessage(fit)))
      } else {
        saveRDS(list(status = "unexpected", model = code, condition_id = cid, rep_id = rid), fit_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – unexpected", code, cid, rid))
      }
    }

    msgs
  }

  ## -- Run with mclapply (forking) -------------------------------------------
  options(mc.cores = 1) # Each Stan fit uses 1 core internally

  results <- parallel::mclapply(
    paths,
    fit_one_dataset,
    mc.cores = cores_outer,
    mc.preschedule = FALSE # Better load balancing for varying runtimes
  )

  ## -- write log -------------------------------------------------------------
  log_entries <- unlist(lapply(results, function(x) {
    if (inherits(x, "try-error")) {
      paste("MCLAPPLY ERROR:", x)
    } else if (is.character(x)) {
      x
    } else {
      NULL
    }
  }))

  if (length(log_entries)) {
    con <- file(log_file, "a")
    writeLines(log_entries, con)
    close(con)
  }

  # Count results
  n_ok <- sum(grepl(": OK", log_entries))
  n_skip <- sum(grepl(": skip", log_entries))
  n_fail <- sum(grepl(": FAIL", log_entries))

  message("\nFinished processing ", length(paths), " datasets.")
  message("  OK: ", n_ok, "  Skip: ", n_skip, "  Fail: ", n_fail)
  message("Log: ", log_file)
}
