###########################################################################
# fit_models.R  –  Fit Skew‑Normal + Normal copula VAR(1) models
# -------------------------------------------------------------------------
# • Parallel over DATASETS via foreach / doParallel  (cores_outer)
# • Within each dataset, Stan chain parallelism controlled by
#   options(mc.cores) ⇢ set with `parallel_chains`
###########################################################################

library(rstan)
library(dplyr)
library(stringr)
library(doParallel)
library(foreach)

fit_models_ml <- function(data_dir,
                          fits_dir,
                          stan_dir,
                          chains = 4,
                          iter = 2000,
                          warmup = 1000,
                          adapt_delta = 0.9,
                          max_treedepth = 12,
                          parallel_chains = TRUE,
                          cores_outer = max(1, parallel::detectCores() - 2)) {
  dir.create(fits_dir, FALSE, TRUE)

  # ---------- cached compilation ----------------------------------------
  cached_model <- function(name, stan_file) {
    rds <- file.path(stan_dir, paste0(name, ".rds"))
    if (file.exists(rds)) {
      readRDS(rds)
    } else {
      mod <- stan_model(stan_file, model_name = name)
      saveRDS(mod, rds)
      mod
    }
  }
  mods <- list(
    SNG = cached_model(
      "SkewN_Gauss_ML",
      file.path(stan_dir, "model_SNG_ml.stan")
    ),
    NG = cached_model(
      "Normal_Gauss_ML",
      file.path(stan_dir, "model_NG_ml.stan")
    )
  )

  # ---------- Stan mc.cores (per dataset) --------------------------------
  options(mc.cores = if (parallel_chains) chains else 1)

  # ---------- helper: convert list‑data to y[N,T,2] ----------------------
  make_y <- function(obj) {
    N <- obj$N
    T <- obj$T
    y <- array(NA_real_, dim = c(N, T, 2))
    for (i in seq_len(N)) {
      tmp <- obj$data |>
        filter(i == !!i) |>
        arrange(t)
      y[i, , 1] <- tmp$y1
      y[i, , 2] <- tmp$y2
    }
    list(N = N, T = T, y = y)
  }

  # ---------- worker for one dataset ------------------------------------
  fit_dataset <- function(dat_path) {
    ids <- str_match(basename(dat_path), "cond(\\d+)_rep(\\d+)")
    cid <- as.integer(ids[2])
    rid <- as.integer(ids[3])

    obj <- readRDS(dat_path)
    sdata <- make_y(obj)

    log_vec <- character()
    for (code in names(mods)) {
      out_file <- file.path(
        fits_dir,
        sprintf("fit_%s_cond%03d_rep%03d.rds", code, cid, rid)
      )
      if (file.exists(out_file)) {
        log_vec <- c(log_vec, paste0(code, "=skip"))
        next
      }

      fit <- sampling(mods[[code]],
        data = sdata,
        chains = chains, iter = iter, warmup = warmup,
        seed = 1000 + cid * 10 + rid,
        control = list(
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth
        ),
        refresh = 0
      )
      saveRDS(fit, out_file)
      log_vec <- c(log_vec, paste0(code, "=OK"))
    }
    paste(basename(dat_path), ":", paste(log_vec, collapse = "; "))
  }

  # ---------- launch parallel fits --------------------------------------
  data_files <- list.files(data_dir,
    "^sim_data_ml_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(data_files)) {
    stop("No simulation files found in ", data_dir)
  }

  cl <- makeCluster(cores_outer)
  registerDoParallel(cl)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
  })

  message(
    sprintf(
      "Fitting %d datasets on %d outer core(s); ",
      length(data_files), cores_outer
    ),
    appendLF = FALSE
  )
  if (parallel_chains) {
    message(sprintf("%d chains per dataset (mc.cores=%d).", chains, chains))
  } else {
    message(sprintf("%d chains per dataset (mc.cores=1).", chains))
  }

  logs <- foreach(
    f = data_files,
    .packages = c("rstan", "dplyr", "stringr"),
    .errorhandling = "pass"
  ) %dopar% fit_dataset(f)

  cat(logs, sep = "\n")
  message(">>> Stan fitting complete.")
}
