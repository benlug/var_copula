##########################################################
# fit_models.R
#
# Fits VAR(1) copula models using single-level Stan code.
# For each simulated dataset, fits:
#   1. The correctly specified model (based on DGP info).
#   2. The standard Normal-Gaussian model (as a baseline).
# Simplified naming and reduced verbosity.
##########################################################

# --- Libraries ---
library(rstan)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(stringr)

# --- Null default operator ---
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Helper Function for Fitting WORKER ---
# log_file: optional path to append Stan messages and errors
fit_specific_models_worker <- function(data_filepath,
                                       model_definitions,
                                       fits_dir,
                                       stan_iter, stan_warmup, stan_chains,
                                       stan_adapt_delta, stan_max_treedepth,
                                       debug_mode = FALSE,
                                       log_file = NULL) {
  # --- Packages needed ---
  require(rstan)
  require(dplyr)
  require(tidyr)
  require(stats)
  require(stringr)

  # --- Parse filename ---
  fname_base <- basename(data_filepath)
  match_res <- str_match(fname_base, "sim_data_cond(\\d+)_rep(\\d+)\\.rds") # Simplified pattern
  if (is.na(match_res[1, 1])) {
    return(paste("ERROR: Cannot parse", fname_base))
  }
  cond_num <- as.integer(match_res[1, 2])
  rep_num <- as.integer(match_res[1, 3])

  # --- Load Data & Get DGP Info ---
  sim_dat <- tryCatch(readRDS(data_filepath), error = function(e) {
    NULL
  })
  if (is.null(sim_dat) || !all(c("T", "data", "correct_fit_model_code") %in% names(sim_dat))) {
    return(paste("ERROR: Invalid data structure in", fname_base))
  }
  T_val <- sim_dat$T
  df <- sim_dat$data
  correct_model_code <- sim_dat$correct_fit_model_code
  standard_model_code <- "NG"
  models_to_fit_codes <- unique(c(correct_model_code, standard_model_code))


  # --- Prepare Stan Data ---
  y_stan_list <- vector("list", T_val)
  tryCatch(
    {
      df_sorted <- df %>% dplyr::arrange(t)
      if (nrow(df_sorted) != T_val || anyNA(df_sorted$y1) || anyNA(df_sorted$y2) || !all(is.finite(c(df_sorted$y1, df_sorted$y2)))) stop("Data invalid")
      for (t_loop in 1:T_val) {
        y_stan_list[[t_loop]] <- c(df_sorted$y1[t_loop], df_sorted$y2[t_loop])
      }
    },
    error = function(e) {
      return(paste("ERROR: Data prep failed", fname_base, ":", conditionMessage(e)))
    }
  )
  # data list for stan
  stan_data <- list(T = T_val, y = y_stan_list)

  # --- random starting values for stan ---
  generate_inits_sl <- function(model_code) {
    base_list <- list(
      mu = stats::rnorm(2, 0, 0.1), phi11 = stats::runif(1, -0.5, 0.5), phi12 = stats::runif(1, -0.3, 0.3),
      phi21 = stats::runif(1, -0.3, 0.3), phi22 = stats::runif(1, -0.5, 0.5)
    )
    sd_y1 <- max(0.1, stats::sd(df$y1, na.rm = TRUE))
    sd_y2 <- max(0.1, stats::sd(df$y2, na.rm = TRUE))
    if (grepl("^N", model_code)) { # Normal Margins
      base_list$sigma <- c(max(0.01, sd_y1 * stats::runif(1, 0.5, 1.5)), max(0.01, sd_y2 * stats::runif(1, 0.5, 1.5)))
    } else { # Skew-Normal Margins
      base_list$omega <- c(max(0.01, sd_y1 * stats::runif(1, 0.5, 1.5)), max(0.01, sd_y2 * stats::runif(1, 0.5, 1.5)))
      base_list$alpha <- stats::rnorm(2, 0, 1.0)
    }
    if (grepl("G$", model_code)) {
      base_list$rho <- stats::runif(1, -0.5, 0.5) # Gaussian Copula
    } else {
      base_list$theta <- abs(stats::rnorm(1, 1, 0.5)) + 0.1
    } # Clayton Copula
    if (any(!sapply(base_list, function(x) all(is.finite(x))))) stop("Non-finite inits")
    return(base_list)
  }

  # --- Fit Required Models ---
  stan_control <- list(adapt_delta = stan_adapt_delta, max_treedepth = stan_max_treedepth)
  common_seed_base <- 2025 * cond_num + rep_num
  fit_results_summary <- list()

  for (code in models_to_fit_codes) {
    model_info <- model_definitions[[code]]
    if (is.null(model_info)) {
      fit_results_summary[[code]] <- "Model def missing"
      next
    }
    compiled_model <- model_info$compiled_model
    model_name_short <- model_info$short_name
    out_file <- file.path(fits_dir, sprintf("fit_%s_cond%03d_rep%03d.rds", model_name_short, cond_num, rep_num)) # Simplified name

    if (file.exists(out_file)) {
      fit_results_summary[[code]] <- paste(code, "skip")
      next
    }

    init_list <- tryCatch(
      {
        lapply(1:stan_chains, function(id) generate_inits_sl(code))
      },
      error = function(e) {
        list(msg = paste("Init Err:", conditionMessage(e)))
      }
    )
    if (!is.list(init_list) || !is.null(init_list$msg)) {
      fit_results_summary[[code]] <- paste(code, "FAIL(Init)")
      next
    }

    fit_obj <- NULL
    fit_error_msg <- NULL
    sampling_msgs <- NULL
    start_time <- Sys.time()
    tryCatch(
      {
        current_refresh <- if (debug_mode) 500 else 0
        # run stan to fit the model; capture warnings/messages (e.g. divergences)
        sampling_msgs <- capture.output(
          fit_obj <- rstan::sampling(
            object = compiled_model, data = stan_data, chains = stan_chains,
            iter = stan_iter, warmup = stan_warmup,
            seed = common_seed_base + which(models_to_fit_codes == code),
            init = init_list, refresh = current_refresh, control = stan_control
          ),
          type = "message"
        )
        saveRDS(fit_obj, file = out_file)
      },
      error = function(e) {
        fit_error_msg <<- paste("Stan Err:", conditionMessage(e))
      }
    )
    if (!is.null(log_file)) {
      # Save Stan warnings (e.g., divergent transitions) alongside errors
      log_conn <- file(log_file, "a")
      if (length(sampling_msgs) > 0) writeLines(sampling_msgs, log_conn)
      if (!is.null(fit_error_msg)) writeLines(fit_error_msg, log_conn)
      close(log_conn)
    }
    end_time <- Sys.time()
    elapsed_time <- difftime(end_time, start_time, units = "mins")

    if (!is.null(fit_error_msg)) {
      fit_results_summary[[code]] <- paste(code, "FAIL(", fit_error_msg, ")")
    } else {
      samples_ok <- FALSE
      if (!is.null(fit_obj)) {
        summary_res <- tryCatch(summary(fit_obj)$summary, error = function(e) NULL)
        if (!is.null(summary_res) && nrow(summary_res) > 0 && all(is.finite(summary_res[, "mean"]))) samples_ok <- TRUE
      }
      fit_results_summary[[code]] <- if (samples_ok) sprintf("%s OK(%.1fm)", code, elapsed_time) else paste(code, "FAIL(Samples)")
    }
  } # End loop over models

  return(paste("Processed:", fname_base, "-", paste(fit_results_summary, collapse = "; ")))
}

# --- Main Fitting Function ---
fit_var1_copula_models <- function(
    data_dir, fits_dir, stan_models_dir, sim_conditions_file,
    stan_iter = 4000, stan_warmup = 2000, stan_chains = 4,
    stan_adapt_delta = 0.99, stan_max_treedepth = 12,
    num_cores = parallel::detectCores() %||% 1,
    debug_mode = FALSE, debug_file_limit = NULL,
    log_file = file.path(fits_dir, "stan_fit_errors.log")) {
  # log_file: path to save Stan fitting errors and warnings; overwritten each run
  # Stan Options & Core Management
  rstan::rstan_options(auto_write = TRUE)
  original_mc_cores <- getOption("mc.cores", 1L)
  cores_to_set <- 1
  if (debug_mode) {
    cores_for_chains <- min(stan_chains, parallel::detectCores(logical = FALSE) %||% 1L)
    cores_to_set <- max(1L, cores_for_chains)
    num_cores <- 1
  }
  cat(sprintf("--- using %d core(s) ---\n", cores_to_set))
  options(mc.cores = cores_to_set)
  on.exit(options(mc.cores = original_mc_cores), add = TRUE)

  # Validate Dirs & Files
  if (!dir.exists(data_dir)) stop("Data directory not found: ", data_dir)
  if (!dir.exists(stan_models_dir)) stop("Stan models directory not found: ", stan_models_dir)
  if (!dir.exists(fits_dir)) dir.create(fits_dir, recursive = TRUE)
  if (!file.exists(sim_conditions_file)) stop("Sim conditions file not found: ", sim_conditions_file)

  if (!is.null(log_file)) {
    try(suppressWarnings(file.remove(log_file)), silent = TRUE)
  }

  # Stan Model Files
  model_files <- list(
    NG = file.path(stan_models_dir, "model_NG_sl.stan"), NC = file.path(stan_models_dir, "model_NC_sl.stan"),
    SG = file.path(stan_models_dir, "model_SNG_sl.stan"), SC = file.path(stan_models_dir, "model_SNC_sl.stan")
  )
  missing_files <- model_files[!sapply(model_files, file.exists)]
  if (length(missing_files) > 0) stop("Stan file(s) not found: ", paste(missing_files, collapse = ", "))

  # compile each stan model once
  cat("compiling stan models...\n")
  compiled_models_list <- list()
  tryCatch(
    {
      compiled_models_list[["NG"]] <- list(short_name = "NG", compiled_model = stan_model(model_files$NG, model_name = "Norm_Gauss"))
      compiled_models_list[["NC"]] <- list(short_name = "NC", compiled_model = stan_model(model_files$NC, model_name = "Norm_Clayton"))
      compiled_models_list[["SG"]] <- list(short_name = "SG", compiled_model = stan_model(model_files$SG, model_name = "SN_Gauss"))
      compiled_models_list[["SC"]] <- list(short_name = "SC", compiled_model = stan_model(model_files$SC, model_name = "SN_Clayton"))
    },
    error = function(e) {
      msg <- paste("Stan compile failed:", e$message)
      if (!is.null(log_file)) try(write(msg, file = log_file, append = TRUE), silent = TRUE)
      stop(msg)
    }
  )
  cat("Stan models compiled.\n")

  # Get Data Files
  all_data_files <- list.files(data_dir, pattern = "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE) # Simplified pattern
  if (length(all_data_files) == 0) {
    cat("No sim_data files found in: ", data_dir, "\n")
    return(invisible(NULL))
  }

  # Apply Debug Limit
  data_files_to_process <- all_data_files
  if (debug_mode && !is.null(debug_file_limit) && debug_file_limit > 0) {
    data_files_to_process <- all_data_files[1:min(debug_file_limit, length(all_data_files))]
  }
  cat(sprintf("found %d simulation files\n", length(data_files_to_process)))
  if (length(data_files_to_process) == 0) {
    cat("No data files selected.\n")
    return(invisible(NULL))
  }

  # Setup Parallel Backend
  cl <- NULL
  run_engine <- NULL
  n_cores_actual <- 1
  if (!debug_mode) {
    n_cores_actual <- max(1, num_cores)
    cat(sprintf("starting parallel cluster with %d cores\n", n_cores_actual))
    if (!requireNamespace("doParallel", quietly = T) || !requireNamespace("foreach", quietly = T)) stop("Packages 'doParallel'/'foreach' needed.")
    if (foreach::getDoParRegistered()) try(doParallel::stopImplicitCluster(), silent = TRUE)
    cl <- tryCatch(parallel::makeCluster(n_cores_actual), error = function(e) stop("Cluster error: ", e$message))
    doParallel::registerDoParallel(cl)
    on.exit(
      {
        if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
        try(doParallel::registerDoSEQ(), silent = TRUE)
      },
      add = TRUE
    )
    run_engine <- foreach::`%dopar%`
  } else {
    run_engine <- foreach::`%do%`
  }

  # Run Loop
  cat(sprintf("starting stan sampling (%s)\n", if (debug_mode) "seq files" else "parallel files"))
  start_time_loop <- Sys.time()

  if (!exists("fit_specific_models_worker", mode = "function")) stop("Worker function not found.")
  vars_to_export <- c(
    "compiled_models_list", "fits_dir", "stan_iter", "stan_warmup", "stan_chains",
    "stan_adapt_delta", "stan_max_treedepth", "debug_mode", "fit_specific_models_worker",
    "log_file"
  )

  foreach_obj <- foreach::foreach(
    f = data_files_to_process, .packages = c("rstan", "dplyr", "tidyr", "stats", "stringr"),
    .export = vars_to_export, .errorhandling = "pass", .verbose = FALSE
  )

  progress_bar <- NULL
  progress_index <- 0
  if (debug_mode) {
    progress_bar <- utils::txtProgressBar(min = 0, max = length(data_files_to_process), style = 3)
  }

  results_list <- run_engine(foreach_obj, {
    res <- fit_specific_models_worker(
      data_filepath = f, model_definitions = compiled_models_list, fits_dir = fits_dir,
      stan_iter = stan_iter, stan_warmup = stan_warmup, stan_chains = stan_chains,
      stan_adapt_delta = stan_adapt_delta, stan_max_treedepth = stan_max_treedepth,
      debug_mode = debug_mode, log_file = log_file
    )
    if (!is.null(progress_bar)) {
      progress_index <<- progress_index + 1
      utils::setTxtProgressBar(progress_bar, progress_index)
    }
    res
  })
  if (!is.null(progress_bar)) close(progress_bar)
  end_time_loop <- Sys.time()
  cat(sprintf("\nLoop finished. Time: %.2f mins.\n", difftime(end_time_loop, start_time_loop, units = "mins")))

  # Process Results Summary
  errors <- list()
  status_messages <- character(length(results_list))
  for (i in seq_along(results_list)) {
    res <- results_list[[i]]
    current_file <- data_files_to_process[i]
    if (inherits(res, "error")) {
      err_msg <- paste("FATAL Worker Error:", conditionMessage(res))
      errors[[length(errors) + 1]] <- list(file = basename(current_file), error = err_msg)
      status_messages[i] <- err_msg
    } else if (is.character(res)) {
      status_messages[i] <- res
      if (grepl("ERROR:|FAILED", res)) {
        is_fatal <- any(sapply(errors, function(e) e$file == basename(current_file) & grepl("FATAL", e$error)))
        if (!is_fatal) errors[[length(errors) + 1]] <- list(file = basename(current_file), error = res)
      }
    } else {
      err_msg <- paste("Unknown result type:", basename(current_file))
      errors[[length(errors) + 1]] <- list(file = basename(current_file), error = err_msg)
      status_messages[i] <- err_msg
    }
  }

  cat("--- fitting summary ---\n")
  cat(status_messages, sep = "\n")
  total_processed <- length(data_files_to_process)
  error_count <- length(errors)
  cat(sprintf("Summary: %d files attempted, %d reported errors.\n", total_processed, error_count))
  if (error_count > 0) {
    cat("Error Details (Max 50):\n")
    error_summary <- list()
    for (err in errors) {
      fname <- err$file
      current_msg <- err$error
      severity <- ifelse(grepl("FATAL", current_msg), 3, ifelse(grepl("ERROR|FAILED", current_msg), 2, 1))
      if (!fname %in% names(error_summary) || severity > error_summary[[fname]]$severity) {
        error_summary[[fname]] <- list(msg = current_msg, severity = severity)
      }
    }
    if (length(error_summary) > 0) {
      error_df <- data.frame(
        file = names(error_summary),
        msg = sapply(error_summary, `[[`, "msg"),
        severity = sapply(error_summary, `[[`, "severity"),
        stringsAsFactors = FALSE
      )
      error_df <- error_df[order(-error_df$severity, error_df$file), ]
      printed_count <- 0
      for (i in 1:min(nrow(error_df), 50)) {
        cat(sprintf("  - %s: %s\n", error_df$file[i], error_df$msg[i]))
        printed_count <- i
      }
      if (nrow(error_df) > printed_count) {
        cat(sprintf("  (... %d more errors ...)\n", nrow(error_df) - printed_count))
      }
    }
    if (!is.null(log_file)) {
      log_conn <- file(log_file, open = "a")
      writeLines(sprintf("=== Stan fitting errors (%s) ===", Sys.time()), log_conn)
      for (err in errors) {
        writeLines(sprintf("%s: %s", err$file, err$error), log_conn)
      }
      close(log_conn)
      cat(sprintf("Error log written to %s\n", log_file))
    }
  } else {
    if (!is.null(log_file)) {
      log_conn <- file(log_file, open = "a")
      writeLines(sprintf("=== Stan fitting run (%s): no errors ===", Sys.time()), log_conn)
      close(log_conn)
    }
  }

  cat("\n--- Fitting script completed. ---\n")
  invisible(NULL)
}
