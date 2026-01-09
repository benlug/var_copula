#!/usr/bin/env Rscript
###########################################################################
# fit_models.R – Study 2 (Gamma extension)
#
# Fits the two comparison models:
#   • NG : Normal margins + Gaussian copula
#   • GG : Gamma margins (shape estimated) + Gaussian copula
#
# Key features
#   • Robust skipping/resume: only skips if an existing .rds is a valid fit
#   • Model cache keyed by Stan-source hash (auto recompile on change)
#   • Parallel outer loop over datasets; chains run sequentially inside workers
#   • Optional lightweight fit artifacts (posterior summaries) for simulation work
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

  save_level <- match.arg(save_level)

  # Parameters needed downstream (analysis/visualization)
  CORE_PARS <- c("mu[1]", "mu[2]", "phi11", "phi12", "phi21", "phi22", "rho")
  WANT_FULL <- list(
    NG = c(CORE_PARS, "sigma[1]", "sigma[2]"),
    GG = c(CORE_PARS, "sigma_gam[1]", "sigma_gam[2]", "shape_gam")
  )

  # Injective seed mapping
  seed_key <- function(cid, rid, chain = 0L) {
    as.integer(cid) * 100000L + as.integer(rid) * 10L + as.integer(chain)
  }

  is_valid_fit_file <- function(x) {
    if (inherits(x, "stanfit")) {
      return(length(x@sim) > 0 && !is.null(x@sim$iter) && x@sim$iter > 0 &&
        length(x@sim$samples) > 0)
    }
    if (is.list(x) && identical(x$type, "stan_fit_summary_v1")) {
      return(is.data.frame(x$summary) && nrow(x$summary) > 0)
    }
    FALSE
  }

  extract_fit_summary <- function(fit, want_full) {
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

  ## -- compile/load cached models (hash by Stan source) ------------------
  model_cache <- function(path, tag) {
    if (!file.exists(path)) return(NULL)

    code <- readChar(path, file.info(path)$size)
    h <- digest::digest(code, algo = "xxhash64")
    rds <- file.path(stan_dir, sprintf("%s_%s.rds", tag, substr(h, 1, 8)))

    # prune older cached versions for this tag
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
    GG = model_cache(file.path(stan_dir, "model_GG_sl.stan"), "Gamma_Gauss")
  )
  models <- Filter(Negate(is.null), models)
  if (!length(models)) stop("No models loaded successfully. Check Stan files.")
  message("Successfully loaded models: ", paste(names(models), collapse = ", "))

  ## -- init helpers ------------------------------------------------------
  make_init_NG <- function() {
    list(
      mu = rnorm(2, 0, 0.1),
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      sigma = runif(2, 0.8, 1.2),
      rho = runif(1, -0.5, 0.5)
    )
  }

  # Data-driven init for GG (helps feasibility + speeds adaptation)
  make_init_GG <- function(y_data, chain_id, skew_dir, shape_true = NULL) {
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
        shape_gam = if (is.finite(shape_true) && shape_true > 0) shape_true else 2,
        rho = runif(1, -0.5, 0.5)
      )
    }

    if (inherits(f1, "try-error") || inherits(f2, "try-error")) return(fallback())
    if (!all(is.finite(coef(f1))) || !all(is.finite(coef(f2)))) return(fallback())

    res <- cbind(residuals(f1), residuals(f2))

    # Use true shape (from simulation) if available; otherwise a conservative default.
    k0 <- if (is.finite(shape_true) && shape_true > 0) shape_true else 2
    # Small jitter to avoid identical inits across chains.
    shape_init <- exp(rnorm(1, log(k0), 0.05))

    list(
      mu    = c(coef(f1)[1], coef(f2)[1]),
      phi11 = coef(f1)["y1_lag"] %||% 0.3,
      phi12 = coef(f1)["y2_lag"] %||% 0.1,
      phi21 = coef(f2)["y1_lag"] %||% 0.1,
      phi22 = coef(f2)["y2_lag"] %||% 0.3,
      eta   = rep(log(0.2), 2),
      shape_gam = shape_init,
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

    # Skew directions (mirroring) encoded in the simulated true parameters
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    m1_mirror <- ds$true_params$margin1$mirror %||% FALSE
    m2_mirror <- ds$true_params$margin2$mirror %||% FALSE
    skew_directions <- c(if (m1_mirror) -1 else 1, if (m2_mirror) -1 else 1)

    # True Gamma shape from DGP (for init only)
    shape_true <- ds$true_params$margin1$shape %||% NA_real_

    msgs <- character()

    for (code in c("NG", "GG")) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))
      want_full <- WANT_FULL[[code]] %||% CORE_PARS

      # Skip if valid
      if (file.exists(fit_path)) {
        chk <- try(readRDS(fit_path), silent = TRUE)
        if (is_valid_fit_file(chk)) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
          next
        } else {
          try(unlink(fit_path), silent = TRUE)
        }
      }

      sdat <- sdat_base
      if (code == "GG") {
        sdat$skew_direction <- skew_directions
        inits <- lapply(seq_len(chains), function(k) {
          make_init_GG(y_matrix, seed_key(cid, rid, k), skew_directions, shape_true = shape_true)
        })
      } else {
        set.seed(SEED_BASE_R + seed_key(cid, rid, 0L))
        inits <- replicate(chains, make_init_NG(), simplify = FALSE)
      }

      seeds_stan <- SEED_BASE_STAN + seed_key(cid, rid, 0L) + 0:(chains - 1)

      if (isTRUE(save_init_sidecar)) {
        try(saveRDS(
          list(
            meta = list(model = code, condition_id = cid, rep_id = rid, time = Sys.time()),
            seeds_stan = seeds_stan,
            init = inits
          ),
          file.path(side_dir, sprintf("INIT_%s_cond%03d_rep%03d.rds", code, cid, rid))
        ), silent = TRUE)
      }

      warn_msgs <- character()
      fit <- withCallingHandlers(
        tryCatch(
          rstan::sampling(
            models[[code]],
            data = sdat,
            chains = chains,
            iter = iter,
            warmup = warmup,
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

        if (save_level == "stanfit") {
          saveRDS(fit, fit_path, compress = fit_compress)
        } else {
          sm_df <- extract_fit_summary(fit, want_full)
          if (is.null(sm_df) || !nrow(sm_df)) {
            # fallback: save full fit so we don't lose information
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
        fail_path <- file.path(fits_dir, sprintf("fit_MISC_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail_path)
        msgs <- c(msgs, sprintf("%s cond%03d rep%03d : FAIL – unexpected return type", code, cid, rid))
      }
    }

    msgs
  }

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
    on.exit({
      if (inherits(con, "connection")) {
        io <- try(isOpen(con), silent = TRUE)
        if (!inherits(io, "try-error") && io) close(con)
      }
    }, add = TRUE)
    try(writeLines(log_entries, con), silent = TRUE)
  }

  message(
    "\nFinished processing ", length(paths), " data sets.",
    "\nLog      → ", log_file,
    "\nSide-cars → ", side_dir
  )
  if (errors_occurred) {
    message("!!! Errors occurred during parallel execution. Check the log file.")
  }
}
