# fit_models_ml.R — fits EG_MLmin & NG_MLmin to ML data
fit_var1_copula_models_ml <- function(data_dir, fits_dir, stan_dir, results_dir,
                                      chains = 4, iter = 4000, warmup = 2000,
                                      adapt_delta = 0.997, max_treedepth = 15,
                                      cores_outer = 2,
                                      start_condition = 1, start_rep = 1) {
  suppressPackageStartupMessages({
    library(rstan)
    library(stringr)
    library(digest)
    library(doParallel)
    library(foreach)
    library(dplyr)
  })

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  dir.create(fits_dir, FALSE, TRUE)
  side_dir <- file.path(results_dir, "init_sidecar_ml")
  dir.create(side_dir, FALSE, TRUE)
  log_file <- file.path(fits_dir, "stan_fit_ml.log")

  model_cache <- function(path, tag) {
    if (!file.exists(path)) {
      return(NULL)
    }
    code <- readChar(path, file.info(path)$size)
    h <- digest(code, algo = "xxhash64")
    rds <- file.path(stan_dir, sprintf("%s_%s.rds", tag, substr(h, 1, 8)))
    old <- list.files(stan_dir, paste0("^", tag, "_[0-9a-f]{8}\\.rds$"), full.names = TRUE)
    old <- setdiff(old, rds)
    if (length(old)) try(unlink(old), silent = TRUE)
    if (file.exists(rds)) {
      message("Loading ", tag, " [", substr(h, 1, 8), "]")
      readRDS(rds)
    } else {
      message("Compiling ", tag, " [", substr(h, 1, 8), "]")
      m <- rstan::stan_model(path, model_name = paste0(tag, "_", substr(h, 1, 8)))
      saveRDS(m, rds)
      m
    }
  }

  models <- list(
    EG_MLmin = model_cache(file.path(stan_dir, "model_EG_ml_min.stan"), "Exp_Gauss_MLmin"),
    NG_MLmin = model_cache(file.path(stan_dir, "model_NG_ml_min.stan"), "Normal_Gauss_MLmin")
  )
  models <- Filter(Negate(is.null), models)
  if (!length(models)) stop("No ML models compiled.")

  # enumerate datasets
  paths <- list.files(data_dir, "^sim_dataML_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(paths)) {
    message("No ML data files found.")
    return(invisible())
  }
  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  keep <- (as.integer(meta[, 2]) > start_condition) |
    (as.integer(meta[, 2]) == start_condition & as.integer(meta[, 3]) >= start_rep)
  paths <- paths[keep]
  if (!length(paths)) {
    message("Nothing to fit (requested cells done).")
    return(invisible())
  }

  # parallel cluster
  if (is.null(cores_outer) || !is.numeric(cores_outer) || cores_outer < 1) cores_outer <- 1
  message("Starting parallel cluster with ", cores_outer, " workers.")
  cl <- parallel::makeCluster(cores_outer)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  })

  options(mc.cores = 1)

  log_vec <- foreach::foreach(
    p = paths, .packages = c("rstan", "digest", "stringr", "dplyr"),
    .errorhandling = "pass"
  ) %dopar% {
    m <- stringr::str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(m[2])
    rid <- as.integer(m[3])

    ds <- try(readRDS(p), silent = TRUE)
    if (inherits(ds, "try-error")) {
      return(sprintf("ERROR cond%03d rep%03d : cannot read.", cid, rid))
    }

    N <- ds$N
    Tn <- ds$T
    sp <- split(ds$data, ds$data$i)
    if (length(sp) != N) {
      return(sprintf("ERROR cond%03d rep%03d : bad split.", cid, rid))
    }

    y_arr <- array(NA_real_, dim = c(N, Tn, 2))
    for (ii in seq_len(N)) {
      y_arr[ii, , 1] <- sp[[ii]]$y1
      y_arr[ii, , 2] <- sp[[ii]]$y2
    }
    sdat_base <- list(N = N, T = Tn, y = y_arr)

    m1_mirror <- ds$true_params$margin1$mirror %||% FALSE
    m2_mirror <- ds$true_params$margin2$mirror %||% FALSE

    .flag <- function(x) isTRUE(as.logical(x)[1])

    skew_dir <- c(
      if (.flag(m1_mirror)) -1 else 1,
      if (.flag(m2_mirror)) -1 else 1
    )

    msgs <- character()

    for (code in names(models)) {
      fit_path <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid))
      # skip if a valid stanfit exists
      if (file.exists(fit_path)) {
        chk <- try(readRDS(fit_path), silent = TRUE)
        if (inherits(chk, "stanfit") && length(chk@sim) > 0 && length(chk@sim$samples) > 0) {
          msgs <- c(msgs, sprintf("%s cond%03d rep%03d : skip (exists)", code, cid, rid))
          next
        } else {
          try(unlink(fit_path), silent = TRUE)
        }
      }

      sdat <- sdat_base
      if (code == "EG_MLmin") sdat$skew_direction <- skew_dir

      # inits: robust, simple
      make_init <- function(code) {
        base <- list(
          mu_bar = c(0, 0),
          tau_mu = c(0.2, 0.2),
          z_mu = matrix(0, nrow = N, ncol = 2),
          phi11 = 0.3, phi12 = 0.1, phi21 = 0.1, phi22 = 0.3,
          rho = runif(1, -0.5, 0.5)
        )
        if (code == "EG_MLmin") {
          c(base, list(eta = rep(log(0.2), 2)))
        } else {
          c(base, list(sigma = c(1.0, 1.0)))
        }
      }

      inits <- replicate(chains, make_init(code), simplify = FALSE)
      seeds <- 5e6 + cid * 1000 + rid + 0:(chains - 1)

      fit <- tryCatch(
        rstan::sampling(
          models[[code]],
          data = sdat,
          chains = chains, iter = iter, warmup = warmup,
          seed = seeds[1], init = inits,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
          refresh = 0, cores = 1
        ),
        error = function(e) e
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
      } else {
        fail_path <- file.path(fits_dir, sprintf("fit_ERR_%s_cond%03d_rep%03d.rds", code, cid, rid))
        saveRDS(fit, fail_path)
        msgs <- c(msgs, sprintf(
          "%s cond%03d rep%03d : FAIL – %s", code, cid, rid,
          if (inherits(fit, "error")) conditionMessage(fit) else "unknown"
        ))
      }
    }
    msgs
  }

  # write log
  if (length(log_vec)) {
    con <- file(log_file, "a")
    on.exit(
      {
        try(close(con), silent = TRUE)
      },
      add = TRUE
    )
    for (v in log_vec) if (is.character(v)) try(writeLines(v, con), silent = TRUE)
  }
  message(
    "\nFinished processing ", length(paths), " ML datasets.",
    "\nLog → ", log_file, "\nSide‑cars → ", side_dir
  )
}
