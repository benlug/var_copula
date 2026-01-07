#!/usr/bin/env Rscript
###########################################################################
# fit_models.R  – Robust sampling + caching + EG init (eta param)
# Changes vs prior:
#   • Model cache keyed by .stan source hash (auto recompile on changes)
#   • Do NOT intercept warnings as return values; only catch errors
#   • Skip-if-exists only for valid stanfit; delete stubs and refit
#   • EG init returns 'eta' (log-slack), not 'sigma_exp'
#   • Save error objects under fit_ERR_* (won’t block retries)
#   • Parallel cleanup uses foreach::registerDoSEQ
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
                                   start_rep = 1,
                                   save_level = c("summary", "stanfit"),
                                   fit_compress = "xz",
                                   save_init_sidecar = TRUE) {
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

  # How much to save per fit file:
  #   - "summary" : small RDS with only posterior summaries + key diagnostics (recommended for simulations)
  #   - "stanfit" : full stanfit object (can be very large)
  save_level <- match.arg(save_level)

  # Parameters we actually need for downstream metrics.
  CORE_PARS <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
  WANT_FULL <- list(
    NG = c(CORE_PARS, "sigma[1]", "sigma[2]"),
    SG = c(CORE_PARS, "omega[1]", "omega[2]", "alpha[1]", "alpha[2]"),
    EG = c(CORE_PARS, "sigma_exp[1]", "sigma_exp[2]"),
    GG = c(CORE_PARS, "sigma_gam[1]", "sigma_gam[2]")
  )

  # Injective seed mapping to avoid collisions if the grid expands.
  seed_key <- function(cid, rid, chain = 0L) {
    as.integer(cid) * 100000L + as.integer(rid) * 10L + as.integer(chain)
  }

  is_valid_fit_file <- function(x) {
    if (inherits(x, "stanfit")) {
      return(length(x@sim) > 0 && !is.null(x@sim$iter) && x@sim$iter > 0 && length(x@sim$samples) > 0)
    }
    if (is.list(x) && identical(x$type, "stan_fit_summary_v1")) {
      return(is.data.frame(x$summary) && nrow(x$summary) > 0)
    }
    FALSE
  }

  extract_fit_summary <- function(fit, want_full) {
    # rstan::summary() works on base parameter names; we then filter to exact indexed names.
    base <- unique(sub("\\[.*$", "", want_full))
    keep <- intersect(base, fit@model_pars)
    if (!length(keep)) return(NULL)

    s <- try(rstan::summary(fit, pars = keep)$summary, silent = TRUE)
    if (inherits(s, "try-error")) return(NULL)

    df <- as.data.frame(s, optional = TRUE)
    df$param <- rownames(df)
    rownames(df) <- NULL
    df <- df[df$param %in% want_full, , drop = FALSE]
    df
  }

  side_dir <- file.path(results_dir, "init_sidecar")
  dir.create(fits_dir, FALSE, TRUE)
  if (isTRUE(save_init_sidecar)) dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit.log")

  ## -- compile or load cached models (hash by source) --------------------
  # Compiles to stan_dir/<tag>_<hash8>.rds and prunes old hashes for same tag
  model_cache <- function(path, tag) {
    if (!file.exists(path)) {
      return(NULL)
    }
    code <- readChar(path, file.info(path)$size)
    h <- digest::digest(code, algo = "xxhash64")
    rds <- file.path(stan_dir, sprintf("%s_%s.rds", tag, substr(h, 1, 8)))

    # prune older versions of same tag
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
    SG = model_cache(file.path(stan_dir, "model_SNG_sl.stan"), "SkewN_Gauss"),
    NG = model_cache(file.path(stan_dir, "model_NG_sl.stan"), "Normal_Gauss"),
    EG = model_cache(file.path(stan_dir, "model_EG_sl.stan"), "Exp_Gauss"),
    GG = model_cache(file.path(stan_dir, "model_GG_sl.stan"), "Gamma_Gauss")
  )
  models <- Filter(Negate(is.null), models)
  if (!length(models)) stop("No models loaded successfully. Check Stan files.")
  message("Successfully loaded models: ", paste(names(models), collapse = ", "))

  ## -- helpers: init lists ----------------------------------------------
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
      NULL
    }
  }

  # EG: data-driven init for eta (log slack); uses skew directions
  make_init_EG <- function(y_data, chain_id, skew_dir) {
    set.seed(SEED_BASE_R + chain_id)

    df <- data.frame(y1 = y_data[, 1], y2 = y_data[, 2])
    df$y1_lag <- dplyr::lag(df$y1)
    df$y2_lag <- dplyr::lag(df$y2)

    f1 <- try(lm(y1 ~ y1_lag + y2_lag, data = df), silent = TRUE)
    f2 <- try(lm(y2 ~ y1_lag + y2_lag, data = df), silent = TRUE)

    fallback <- function() {
      list(
        mu = rnorm(2, 0, 0.1),
        phi11 = 0.3, phi12 = 0.1, phi21 = 0.1, phi22 = 0.3,
        eta = rep(log(0.2), 2), # slack ≈ 0.2
        rho = runif(1, -0.5, 0.5)
      )
    }

    if (inherits(f1, "try-error") || inherits(f2, "try-error")) {
      return(fallback())
    }
    if (!all(is.finite(coef(f1))) || !all(is.finite(coef(f2)))) {
      return(fallback())
    }

    res <- cbind(residuals(f1), residuals(f2))
    b <- c(
      max(-skew_dir[1] * res[, 1], na.rm = TRUE),
      max(-skew_dir[2] * res[, 2], na.rm = TRUE)
    )
    list(
      mu    = c(coef(f1)[1], coef(f2)[1]),
      phi11 = coef(f1)["y1_lag"] %||% 0.3,
      phi12 = coef(f1)["y2_lag"] %||% 0.1,
      phi21 = coef(f2)["y1_lag"] %||% 0.1,
      phi22 = coef(f2)["y2_lag"] %||% 0.3,
      eta   = rep(log(0.2), 2), # ETA, not sigma_exp
      rho   = runif(1, -0.5, 0.5)
    )
  }

  # GG: data-driven init for eta (log slack); uses skew directions
  # sigma_gam = b + exp(eta), where b = max(-s * res)/sqrt(shape_gam)
  make_init_GG <- function(y_data, chain_id, skew_dir, shape_gam) {
    set.seed(SEED_BASE_R + chain_id)

    df <- data.frame(y1 = y_data[, 1], y2 = y_data[, 2])
    df$y1_lag <- dplyr::lag(df$y1)
    df$y2_lag <- dplyr::lag(df$y2)

    f1 <- try(lm(y1 ~ y1_lag + y2_lag, data = df), silent = TRUE)
    f2 <- try(lm(y2 ~ y1_lag + y2_lag, data = df), silent = TRUE)

    fallback <- function() {
      list(
        mu = rnorm(2, 0, 0.1),
        phi11 = 0.3, phi12 = 0.1, phi21 = 0.1, phi22 = 0.3,
        eta = rep(log(0.2), 2),
        rho = runif(1, -0.5, 0.5)
      )
    }

    if (inherits(f1, "try-error") || inherits(f2, "try-error")) {
      return(fallback())
    }
    if (!all(is.finite(coef(f1))) || !all(is.finite(coef(f2)))) {
      return(fallback())
    }

    # residuals at the OLS init (only used to compute the feasibility bound b)
    res <- cbind(residuals(f1), residuals(f2))
    sqrt_shape <- sqrt(shape_gam)
    b <- c(
      max(-skew_dir[1] * res[, 1], na.rm = TRUE) / sqrt_shape,
      max(-skew_dir[2] * res[, 2], na.rm = TRUE) / sqrt_shape
    )

    list(
      mu    = c(coef(f1)[1], coef(f2)[1]),
      phi11 = coef(f1)["y1_lag"] %||% 0.3,
      phi12 = coef(f1)["y2_lag"] %||% 0.1,
      phi21 = coef(f2)["y1_lag"] %||% 0.1,
      phi22 = coef(f2)["y2_lag"] %||% 0.3,
      # initialize slightly above feasibility bound via eta (slack = 0.2)
      eta   = rep(log(0.2), 2),
      rho   = runif(1, -0.5, 0.5)
    )
  }

  ## -- enumerate datasets & resume --------------------------------------
  paths <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  if (!length(paths)) {
    message("No data files found in ", data_dir)
    return(invisible())
  }
  if (is.null(meta) || nrow(meta) != length(paths) || ncol(meta) < 3) {
    stop("Error parsing filenames in data directory.")
  }

  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition & as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  if (!length(paths)) {
    message("Nothing to fit (all requested cells done).")
    return(invisible())
  }

  ## -- parallel ----------------------------------------------------------
  if (is.null(cores_outer) || !is.numeric(cores_outer) || cores_outer < 1) cores_outer <- 1
  message("Starting parallel cluster with ", cores_outer, " cores.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    if (requireNamespace("foreach", quietly = TRUE)) try(foreach::registerDoSEQ(), silent = TRUE)
  })

  options(mc.cores = 1) # avoid nested parallelism inside workers

  log_vec <- foreach::foreach(
    p              = paths,
    .packages      = c("rstan", "digest", "stringr", "dplyr"),
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

    # Skew directions from DGP mirrors
    m1_mirror <- ds$true_params$margin1$mirror %||% FALSE
    m2_mirror <- ds$true_params$margin2$mirror %||% FALSE
    skew_directions <- c(if (m1_mirror) -1 else 1, if (m2_mirror) -1 else 1)

    msgs <- character()

    # Decide which models to run
    dgp_type <- ds$true_params$margin1$type %||% "unknown"
    is_exp_dgp <- identical(dgp_type, "exponential")
    is_gam_dgp <- identical(dgp_type, "gamma")

    # Default: run all available models, but only include the correctly-specified
    # positive-support model when it matches the DGP type.
    if (is_exp_dgp) {
      models_to_run <- setdiff(names(models), "GG")
    } else if (is_gam_dgp) {
      models_to_run <- setdiff(names(models), "EG")
    } else {
      models_to_run <- setdiff(names(models), c("EG", "GG"))
    }
    models_to_run <- intersect(models_to_run, names(models))

    for (code in models_to_run) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))

      # Parameters we care about for this model (used for lightweight summaries)
      want_full <- WANT_FULL[[code]] %||% CORE_PARS

      # Skip only if an existing file is a valid completed fit artifact; otherwise refit.
      # If we are in "summary" mode and a legacy full stanfit exists, convert it in-place.
      if (file.exists(fit_path)) {
        chk <- try(readRDS(fit_path), silent = TRUE)
        if (save_level == "summary" && inherits(chk, "stanfit") && is_valid_fit_file(chk)) {
          # Convert existing full stanfit -> lightweight summary (reduces disk immediately, no refit)
          divs_old <- tryCatch(
            {
              sp <- rstan::get_sampler_params(chk, inc_warmup = FALSE)
              if (is.list(sp) && length(sp) > 0) sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L)) else NA_integer_
            },
            error = function(e) NA_integer_
          )

          sm_old <- extract_fit_summary(chk, want_full)
          if (is.data.frame(sm_old) && nrow(sm_old) > 0) {
            fit_light <- list(
              type = "stan_fit_summary_v1",
              model = code,
              condition_id = cid,
              rep_id = rid,
              time_saved = Sys.time(),
              stan_config = list(
                chains = chains, iter = iter, warmup = warmup,
                adapt_delta = adapt_delta, max_treedepth = max_treedepth
              ),
              seeds_stan = NA,
              n_div = divs_old,
              warnings = character(),
              summary = sm_old
            )
            saveRDS(fit_light, fit_path, compress = fit_compress)
            msgs <- c(msgs, sprintf("%s cond%03d rep%03d : CONVERT existing stanfit -> summary", code, cid, rid))
            next
          }
          # If summary extraction fails, fall through to the standard validity check.
        }
        if (is_valid_fit_file(chk)) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
          next
        } else {
          try(unlink(fit_path), silent = TRUE) # remove stub/error so we can refit
        }
      }

      sdat <- sdat_base
      # EG / GG need skew_direction in data and eta init
      if (code %in% c("EG", "GG")) {
        sdat$skew_direction <- skew_directions
        if (code == "GG") {
          # fixed gamma shape passed via data
          shape_gam <- ds$true_params$margin1$shape %||% NA_real_
          if (!is.finite(shape_gam) || shape_gam <= 0) {
            stop("GG model requires ds$true_params$margin1$shape to be a positive number")
          }
          sdat$shape_gam <- shape_gam
          inits <- lapply(seq_len(chains), function(k) {
            make_init_GG(y_matrix, seed_key(cid, rid, k), skew_directions, shape_gam)
          })
        } else {
          # EG
          inits <- lapply(seq_len(chains), function(k) {
            make_init_EG(y_matrix, seed_key(cid, rid, k), skew_directions)
          })
        }
      } else {
        # reproducible inits per dataset
        set.seed(SEED_BASE_R + seed_key(cid, rid, 0L))
        inits <- replicate(chains, make_init(code), simplify = FALSE)
      }

      if (!length(inits) || is.null(inits[[1]])) {
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – init generation failed", code, cid, rid))
        next
      }

      seeds_stan <- SEED_BASE_STAN + seed_key(cid, rid, 0L) + 0:(chains - 1)

      # Save init sidecar (optional; useful for debugging/repro, but can be disk-heavy at scale)
      if (isTRUE(save_init_sidecar)) {
        try(saveRDS(
          list(
            meta = list(model = code, condition_id = cid, rep_id = rid, time = Sys.time()),
            seeds_stan = seeds_stan, init = inits
          ),
          file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
        ), silent = TRUE)
      }

      ## --- sampling (catch errors, NOT warnings) -----------------------
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
        # Divergence count (post-warmup)
        divs <- tryCatch(
          {
            sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
            if (is.list(sp) && length(sp) > 0) sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0L)) else NA_integer_
          },
          error = function(e) NA_integer_
        )

        if (save_level == "stanfit") {
          # Full stanfit (can be very large)
          saveRDS(fit, fit_path, compress = fit_compress)
        } else {
          # Lightweight artifact: posterior summaries for parameters we need + key diagnostics
          want_full <- WANT_FULL[[code]] %||% CORE_PARS
          sm_df <- extract_fit_summary(fit, want_full)
          if (is.null(sm_df) || !nrow(sm_df)) {
            # If summary extraction fails, fall back to saving the full fit (so we don't silently lose information)
            saveRDS(fit, fit_path, compress = fit_compress)
            msgs <- c(msgs, sprintf("%s cond%03d rep%03d : WARN – summary extraction failed; saved full stanfit", code, cid, rid))
          } else {
            fit_light <- list(
              type = "stan_fit_summary_v1",
              model = code,
              condition_id = cid,
              rep_id = rid,
              time_saved = Sys.time(),
              stan_config = list(
                chains = chains, iter = iter, warmup = warmup,
                adapt_delta = adapt_delta, max_treedepth = max_treedepth
              ),
              seeds_stan = seeds_stan,
              n_div = divs,
              warnings = unique(warn_msgs),
              summary = sm_df
            )
            saveRDS(fit_light, fit_path, compress = fit_compress)
          }
        }

        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : OK (%s div)",
          code, cid, rid, ifelse(is.na(divs), "NA", as.character(divs))
        ))
        if (length(warn_msgs)) {
          msgs <- c(msgs, sprintf(
            "%s cond%03d rep%03d : WARN – %s",
            code, cid, rid, paste(unique(warn_msgs), collapse = " | ")
          ))
        }
      } else if (inherits(fit, "error")) {
        fail_path <- file.path(fits_dir, sprintf("fit_ERR_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – %s", code, cid, rid, conditionMessage(fit)))
      } else {
        # Unexpected type (very rare); save separately
        fail_path <- file.path(fits_dir, sprintf("fit_MISC_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – unexpected return type", code, cid, rid))
      }
    } # models loop

    msgs
  } # foreach

  ## -- write log safely --------------------------------------------------
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
  if (errors_occurred) {
    message("!!! Errors occurred during parallel execution. Check the log file.")
  }
}
