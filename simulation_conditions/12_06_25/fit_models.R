###########################################################################
# fit_models.R   (single‑level, only Normal vs. Skew‑Normal)
#
# * Compiles model_NG_sl.stan   and   model_SNG_sl.stan
# * Fits both to every simulated data set, saves    fit_<CODE>_cond###_rep###.rds
# * Stores short run‑log, divergence info etc.
###########################################################################

library(rstan)
library(dplyr)
library(tidyr)
library(stringr)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------
# Worker – fits NG and SG to one .rds data file
# ------------------------------------------------------------------------
fit_specific_models_worker <- function(data_filepath,
                                       model_definitions,
                                       fits_dir,
                                       stan_iter, stan_warmup, stan_chains,
                                       stan_adapt_delta, stan_max_treedepth,
                                       log_file = NULL) {
  fname <- basename(data_filepath)
  m <- str_match(fname, "sim_data_cond(\\d+)_rep(\\d+)\\.rds")
  if (is.na(m[1, 1])) {
    return(paste("ERROR: cannot parse", fname))
  }
  cond_id <- as.integer(m[1, 2])
  rep_id <- as.integer(m[1, 3])

  sim_dat <- tryCatch(readRDS(data_filepath), error = function(e) NULL)
  if (is.null(sim_dat) || !is.list(sim_dat$data)) {
    return(paste("ERROR: bad data in", fname))
  }

  T_val <- sim_dat$T
  y_df <- sim_dat$data %>% arrange(t)
  stan_data <- list(
    T = T_val,
    y = lapply(
      seq_len(T_val),
      function(i) c(y_df$y1[i], y_df$y2[i])
    )
  )

  # models to fit: correct (SG) plus NG
  models_to_fit <- unique(c("SG", "NG"))
  results <- list()

  for (code in models_to_fit) {
    mdl <- model_definitions[[code]]
    outfile <- file.path(
      fits_dir,
      sprintf(
        "fit_%s_cond%03d_rep%03d.rds",
        code, cond_id, rep_id
      )
    )
    if (file.exists(outfile)) {
      results[[code]] <- paste(code, "skip")
      next
    }
    init_fun <- function() {
      list(
        mu = rnorm(2, 0, 0.1),
        phi11 = runif(1, -0.5, 0.5), phi12 = runif(1, -0.3, 0.3),
        phi21 = runif(1, -0.3, 0.3), phi22 = runif(1, -0.5, 0.5),
        sigma = runif(2, 0.5, 1.5),
        rho = runif(1, -0.4, 0.4),
        omega = if (code == "SG") runif(2, 0.5, 1.5) else NULL,
        alpha = if (code == "SG") rnorm(2, 0, 2) else NULL
      )
    }
    t0 <- Sys.time()
    fit <- tryCatch(
      sampling(
        mdl$compiled_model,
        data = stan_data,
        chains = stan_chains,
        iter = stan_iter,
        warmup = stan_warmup,
        seed = 1000 + cond_id * 10 + rep_id,
        control = list(
          adapt_delta = stan_adapt_delta,
          max_treedepth = stan_max_treedepth
        ),
        init = init_fun,
        refresh = 0
      ),
      error = function(e) e
    )
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 2)
    if (inherits(fit, "stanfit")) {
      saveRDS(fit, outfile)
      results[[code]] <- sprintf("%s OK (%.2f min)", code, elapsed)
    } else {
      msg <- paste(code, "FAIL:", fit$message)
      if (!is.null(log_file)) write(msg, file = log_file, append = TRUE)
      results[[code]] <- msg
    }
  }
  paste("Processed", fname, "->", paste(results, collapse = "; "))
}

# ------------------------------------------------------------------------
# Main wrapper
# ------------------------------------------------------------------------
fit_var1_copula_models <- function(data_dir, fits_dir, stan_models_dir,
                                   sim_conditions_file,
                                   stan_iter = 4000, stan_warmup = 2000,
                                   stan_chains = 4,
                                   stan_adapt_delta = 0.9,
                                   stan_max_treedepth = 12,
                                   num_cores = parallel::detectCores() %||% 1,
                                   log_file = file.path(fits_dir, "stan_fit.log")) {
  dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(sim_conditions_file)) {
    stop("conditions file missing")
  }
  print(stan_models_dir)

  # -------------------- compile the two Stan models ----------------------
  model_files <- list(
    NG = file.path(stan_models_dir, "model_NG_sl.stan"),
    SG = file.path(stan_models_dir, "model_SNG_sl.stan")
  )
  model_paths <- unlist(model_files, use.names = FALSE)
  if (!all(file.exists(model_paths))) {
    stop(
      "Stan model file(s) missing: ",
      paste(basename(model_paths)[!file.exists(model_paths)], collapse = ", ")
    )
  }

  cat("Compiling Stan models ...\n")
  compiled_models <- list(
    NG = list(
      short_name = "NG",
      compiled_model = stan_model(model_files$NG,
        model_name = "Normal_Gauss"
      )
    ),
    SG = list(
      short_name = "SG",
      compiled_model = stan_model(model_files$SG,
        model_name = "SkewN_Gauss"
      )
    )
  )

  # -------------------- pick data files ----------------------------------
  data_files <- list.files(data_dir,
    "^sim_data_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(data_files)) stop("no data files found")

  # parallel cluster
  cores <- max(1, num_cores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(
    {
      parallel::stopCluster(cl)
      doParallel::registerDoSEQ()
    },
    add = TRUE
  )

  cat(sprintf(
    "Launching fits (%d files, %d cores)\n",
    length(data_files), cores
  ))

  `%dopar%` <- foreach::`%dopar%`
  res <- foreach::foreach(
    f = data_files,
    .packages = c("rstan", "dplyr", "tidyr", "stringr"),
    .export = c(
      "fit_specific_models_worker",
      "%||%"
    ),
    .errorhandling = "pass"
  ) %dopar% {
    fit_specific_models_worker(
      data_filepath      = f,
      model_definitions  = compiled_models,
      fits_dir           = fits_dir,
      stan_iter          = stan_iter,
      stan_warmup        = stan_warmup,
      stan_chains        = stan_chains,
      stan_adapt_delta   = stan_adapt_delta,
      stan_max_treedepth = stan_max_treedepth,
      log_file           = log_file
    )
  }
  cat(res, sep = "\n")
  cat("All fits finished.\n")
}
