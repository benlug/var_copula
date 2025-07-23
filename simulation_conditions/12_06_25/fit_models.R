###########################################################################
# fit_models.R  – robust Stan fitting (single‑level)
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
      m <- stan_model(path, model_name = tag)
      saveRDS(m, rds)
      m
    }
  }
  mod_SG <- model_cache(
    file.path(stan_dir, "model_SNG_sl.stan"),
    "SkewN_Gauss"
  )
  mod_NG <- model_cache(
    file.path(stan_dir, "model_NG_sl.stan"),
    "Normal_Gauss"
  )
  models <- list(SG = mod_SG, NG = mod_NG)

  ## -- helper: minimal init list per model ------------------------------
  make_init <- function(code) {
    base <- list(
      mu = rnorm(2, 0, 0.1),
      phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
      phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
      sigma = runif(2, 0.5, 2), rho = runif(1, -0.5, 0.5)
    )
    if (code == "SG") {
      c(base,
        omega = runif(2, 0.5, 2),
        alpha = runif(2, -10, 10)
      )
    } else {
      base
    }
  }

  ## -- list data sets & honour resume -----------------------------------
  paths <- list.files(data_dir,
    "^sim_data_cond\\d+_rep\\d+\\.rds$",
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

  options(mc.cores = parallel::detectCores())

  log_vec <- foreach(
    p = paths,
    .packages = c("rstan", "digest", "stringr"),
    .errorhandling = "pass"
  ) %dopar% {
    m <- str_match(basename(p), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(m[2])
    rid <- as.integer(m[3])

    ds <- readRDS(p)
    sdat <- list(
      T = ds$T,
      y = lapply(
        seq_len(ds$T),
        \(i) c(ds$data$y1[i], ds$data$y2[i])
      )
    )

    msgs <- character()

    for (code in names(models)) {
      fit_path <- file.path(
        fits_dir, sprintf(
          "fit_%s_cond%03d_rep%03d.rds",
          code, cid, rid
        )
      )
      if (file.exists(fit_path)) {
        msgs <- c(
          msgs,
          sprintf(
            "%s cond%03d rep%03d : skip (exists)",
            code, cid, rid
          )
        )
        next
      }

      ## robust init --------------
      set.seed(1000 + cid * 10 + rid)
      inits <- replicate(chains, make_init(code), simplify = FALSE)

      sidecar <- list(
        meta = list(
          model = code, condition_id = cid, rep_id = rid,
          time = Sys.time()
        ),
        seeds = 1000 + cid * 10 + rid + 0:(chains - 1),
        init = inits
      )
      saveRDS(
        sidecar,
        file.path(
          side_dir,
          sprintf(
            "INIT_%s_cond%03d_rep%03d.rds",
            code, cid, rid
          )
        )
      )

      ## sampling ------------------
      fit <- tryCatch(
        sampling(models[[code]],
          data = sdat,
          chains = chains, iter = iter, warmup = warmup,
          seed = sidecar$seeds[1],
          init = inits,
          control = list(
            adapt_delta = adapt_delta,
            max_treedepth = max_treedepth
          ),
          refresh = 0
        ),
        error = function(e) e
      )

      saveRDS(fit, fit_path)

      if (inherits(fit, "stanfit")) {
        ndiv <- sum(vapply(
          get_sampler_params(fit, FALSE),
          \(x) sum(x[, "divergent__"]), 0
        ))
        msgs <- c(
          msgs,
          sprintf(
            "%s cond%03d rep%03d : OK (%d div)",
            code, cid, rid, ndiv
          )
        )
      } else {
        msgs <- c(
          msgs,
          sprintf(
            "%s cond%03d rep%03d : FAIL – %s",
            code, cid, rid, conditionMessage(fit)
          )
        )
      }
    }
    msgs
  }

  ## -- write log safely --------------------------------------------------
  log_chr <- unlist(lapply(log_vec, function(x) {
    if (is.character(x)) {
      x
    } else {
      sprintf("WORKER‑ERROR: class %s", class(x)[1])
    }
  }), use.names = FALSE)
  writeLines(log_chr, log_file)
  message(
    "\nFinished ", length(log_chr),
    " model fits. Log → ", log_file,
    "\nSide‑cars → ", side_dir
  )
}
