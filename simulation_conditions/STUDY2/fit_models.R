#!/usr/bin/env Rscript
###########################################################################
# fit_models.R  – (FIXED: Corrected namespace for registerDoSEQ; Robust initialization and cleanup)
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
    # Load necessary libraries and ensure they are available
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
  # Define the helper function within the scope so foreach detects it
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  side_dir <- file.path(results_dir, "init_sidecar")
  dir.create(fits_dir, FALSE, TRUE)
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  ## -- compile or load cached models ------------------------------------
  # The design allows models to be optionally present.
  model_cache <- function(path, tag) {
    rds <- file.path(stan_dir, paste0(tag, ".rds"))
    if (file.exists(path)) {
      if (file.exists(rds)) {
        message("Loading cached model: ", tag)
        return(readRDS(rds))
      } else {
        message("Compiling model: ", tag)
        # Use rstan::stan_model to be explicit
        m <- rstan::stan_model(path, model_name = tag)
        saveRDS(m, rds)
        return(m)
      }
    } else {
      # Silenced: message("Stan file not found for ", tag)
      return(NULL)
    }
  }

  models <- list(
    SG = model_cache(file.path(stan_dir, "model_SNG_sl.stan"), "SkewN_Gauss"),
    NG = model_cache(file.path(stan_dir, "model_NG_sl.stan"), "Normal_Gauss"),
    EG = model_cache(file.path(stan_dir, "model_EG_sl.stan"), "Exp_Gauss")
  )

  # Remove models that didn't load
  models <- Filter(Negate(is.null), models)
  if (length(models) == 0) stop("No models loaded successfully. Check Stan files.")

  # Report which models were successfully loaded
  message("Successfully loaded models: ", paste(names(models), collapse = ", "))


  ## -- helper: minimal init list per model ------------------------------
  make_init <- function(code) {
    base <- list(
      mu = rnorm(2, 0, 0.1),
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      rho = runif(1, -0.5, 0.5)
    )
    if (code == "SG") {
      c(base, list(omega = runif(2, 0.8, 1.5), delta = runif(2, -0.8, 0.8)))
    } else if (code == "NG") {
      c(base, list(sigma = runif(2, 0.8, 1.2)))
    } else {
      return(NULL)
    }
  }

  ## -- helper: INTELLIGENT init for EG model -------------------------------
  # This function provides a "warm start" to avoid support violations.
  make_init_EG <- function(y_data, chain_id) {
    # Use a slightly different seed per chain for random part of init
    set.seed(SEED_BASE_R + chain_id)

    df <- as.data.frame(y_data)
    names(df) <- c("y1", "y2")
    # Use dplyr::lag for safety
    df$y1_lag <- dplyr::lag(df$y1)
    df$y2_lag <- dplyr::lag(df$y2)

    # Fit simple linear models to get plausible starting points
    # Robustness check for lm (e.g., very short T)
    fit1 <- try(lm(y1 ~ y1_lag + y2_lag, data = df), silent = TRUE)
    fit2 <- try(lm(y2 ~ y1_lag + y2_lag, data = df), silent = TRUE)

    # Fallback initialization if lm fails or produces NAs/Infs
    use_fallback <- inherits(fit1, "try-error") || inherits(fit2, "try-error") ||
      !all(is.finite(coef(fit1))) || !all(is.finite(coef(fit2)))

    if (use_fallback) {
      return(list(
        mu = rnorm(2, 0, 0.1),
        phi11 = 0.3, phi12 = 0.1, phi21 = 0.1, phi22 = 0.3,
        sigma_exp = c(1, 1),
        rho = 0.3
      ))
    }

    # Get initial residuals
    res1 <- residuals(fit1)
    res2 <- residuals(fit2)

    # Set sigma_exp to be slightly larger than the max absolute residual
    # This ensures the support condition is met at the start.
    init_sigma1 <- max(abs(res1), na.rm = TRUE) * 1.1 + 0.1
    init_sigma2 <- max(abs(res2), na.rm = TRUE) * 1.1 + 0.1

    # Extract coefficients safely
    coef1 <- coef(fit1)
    coef2 <- coef(fit2)

    list(
      # Initialize mu and phi from the linear model fits
      mu = c(coef1[1], coef2[1]),
      # Use default if NA (can happen with perfect collinearity, though unlikely here)
      phi11 = ifelse(is.na(coef1["y1_lag"]), 0.3, coef1["y1_lag"]),
      phi12 = ifelse(is.na(coef1["y2_lag"]), 0.1, coef1["y2_lag"]),
      phi21 = ifelse(is.na(coef2["y1_lag"]), 0.1, coef2["y1_lag"]),
      phi22 = ifelse(is.na(coef2["y2_lag"]), 0.3, coef2["y2_lag"]),
      # Initialize sigma_exp with a safe margin
      sigma_exp = c(init_sigma1, init_sigma2),
      # Randomly initialize rho
      rho = runif(1, -0.5, 0.5)
    )
  }


  ## -- list data sets & honour resume -----------------------------------
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")

  if (length(paths) == 0) {
    message("No data files found in ", data_dir)
    return(invisible())
  }

  # Ensure meta extraction was successful
  if (is.null(meta) || nrow(meta) != length(paths) || ncol(meta) < 3) {
    stop("Error parsing filenames in data directory.")
  }

  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition &
      as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  if (!length(paths)) {
    message("Nothing to fit (all requested cells done).")
    return(invisible())
  }

  ## -- foreach cluster ---------------------------------------------------
  # Check if cores_outer is valid
  if (is.null(cores_outer) || !is.numeric(cores_outer) || cores_outer < 1) {
    cores_outer <- 1
  }

  message("Starting parallel cluster with ", cores_outer, " cores.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)

  # Robust cleanup: Ensure the cluster is stopped and the backend is reset to sequential
  on.exit({
    # message("Cleaning up parallel environment.") # Optional: message for verbosity
    try(parallel::stopCluster(cl), silent = TRUE)
    # FIX: registerDoSEQ belongs to the 'foreach' package, not 'doParallel'.
    # We call it robustly to ensure the environment is clean after the function exits.
    if (requireNamespace("foreach", quietly = TRUE)) {
      try(foreach::registerDoSEQ(), silent = TRUE)
    }
  })

  # Set internal Stan cores to 1 to avoid nested parallelism
  options(mc.cores = 1)

  # foreach automatically detects the required variables/functions (models, make_init_EG, etc.)
  log_vec <- foreach::foreach(
    p              = paths,
    .packages      = c("rstan", "digest", "stringr", "dplyr"),
    .errorhandling = "pass" # Allows catching errors from workers
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

    # Handle potential missing mirror info gracefully
    m1_mirror <- ds$true_params$margin1$mirror %||% FALSE
    m2_mirror <- ds$true_params$margin2$mirror %||% FALSE
    skew_directions <- c(if (m1_mirror) -1 else 1, if (m2_mirror) -1 else 1)

    msgs <- character()

    # Determine DGP type
    dgp_type <- ds$true_params$margin1$type %||% "unknown"
    is_exponential_dgp <- dgp_type == "exponential"

    # Logic: If exponential DGP, try to run all *loaded* models.
    # If not exponential (e.g., Study 1), run all loaded models except EG.
    models_to_run <- if (is_exponential_dgp) names(models) else setdiff(names(models), "EG")

    # Filter against models that are actually loaded
    models_to_run <- intersect(models_to_run, names(models))

    for (code in models_to_run) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))
      if (file.exists(fit_path)) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
        next
      }

      sdat <- sdat_base

      # Select the correct initialization function
      if (code == "EG") {
        sdat$skew_direction <- skew_directions
        # Use the data-driven init function for the EG model
        inits <- lapply(1:chains, function(chain_id) make_init_EG(y_matrix, chain_id))
      } else {
        # Use the standard random init for NG and SG models
        # Ensure reproducibility within the parallel loop
        set.seed(SEED_BASE_R + cid * 1000 + rid)
        inits <- replicate(chains, make_init(code), simplify = FALSE)
      }

      # Check if initialization generation failed
      if (is.null(inits) || length(inits) == 0 || is.null(inits[[1]])) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – Initialization generation failed", code, cid, rid))
        next
      }

      seeds_stan <- SEED_BASE_STAN + cid * 1000 + rid + 0:(chains - 1)

      # Save initialization sidecar
      try(saveRDS(
        list(
          meta = list(model = code, condition_id = cid, rep_id = rid, time = Sys.time()),
          seeds_stan = seeds_stan, init = inits
        ),
        file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
      ))

      ## --- sampling ----------------------------------------------------
      fit <- tryCatch(
        rstan::sampling(models[[code]],
          data = sdat,
          chains = chains, iter = iter, warmup = warmup,
          seed = seeds_stan[1],
          init = inits,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
          refresh = 0, # Suppress Stan's internal progress messages
          cores = 1 # Ensure single core usage within the worker
        ),
        error = function(e) e,
        warning = function(w) w
      )

      # Save the result (even if it's an error object)
      saveRDS(fit, fit_path)

      if (inherits(fit, "stanfit") && length(fit@sim) > 0) {
        # Check for divergences
        divs <- tryCatch(
          {
            sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
            if (is.list(sp) && length(sp) > 0) {
              sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L))
            } else {
              NA_integer_
            }
          },
          error = function(e) NA_integer_
        )

        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : OK (%d div)", code, cid, rid, divs))
      } else {
        # Handle error or warning results, or empty stanfit objects
        msg <- if (inherits(fit, "error")) conditionMessage(fit) else if (inherits(fit, "warning")) conditionMessage(fit) else "Unknown failure or empty fit"
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – %s", code, cid, rid, msg))
      }
    }
    # Return the log messages for this dataset
    msgs
  }

  ## -- write log safely --------------------------------------------------
  # Process results from foreach loop
  log_entries <- character()
  errors_occurred <- FALSE

  for (result in log_vec) {
    if (inherits(result, "error")) {
      # Handle errors that occurred during the %dopar% execution itself
      log_entries <- c(log_entries, paste("FOREACH ERROR:", conditionMessage(result)))
      errors_occurred <- TRUE
    } else if (is.character(result)) {
      log_entries <- c(log_entries, result)
    }
  }

  if (length(log_entries) > 0) {
    # Append to the log file
    # Use base::file to ensure correct connection handling
    con <- base::file(log_file, "a")
    tryCatch({
      writeLines(log_entries, con)
    }, finally = {
      # Check if connection is valid before attempting to close
      if (inherits(con, "connection")) {
        # Use try(isOpen(con)) for maximum robustness in different R environments
        is_open <- try(isOpen(con), silent = TRUE)
        if (!inherits(is_open, "try-error") && is_open) {
          close(con)
        }
      }
    })
  }

  # The on.exit() handler will execute now, stopping the cluster.

  message(
    "\nFinished processing ", length(paths), " data sets.",
    "\nLog      → ", log_file,
    "\nSide‑cars → ", side_dir
  )
  if (errors_occurred) {
    message("!!! Errors occurred during parallel execution. Check the log file.")
  }
}
