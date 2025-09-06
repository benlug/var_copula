#!/usr/bin/env Rscript
###########################################################################
# inspect_fit.R  –  Quick diagnostics for a single Stan .rds file
#
# Usage (terminal):
#   Rscript inspect_fit.R path/to/fit_file.rds
#
# Usage (interactive):
#   source("inspect_fit.R")
#   inspect_fit("fits/fit_SG_cond001_rep001.rds")
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  # Use posterior package for modern diagnostics (ESS, Rhat)
  library(posterior)
  library(dplyr)
  library(tibble)
})

safe_read <- function(path)
  tryCatch(readRDS(path), error = function(e) e)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
# E-BFMI calculation (Energy Bayesian Fraction of Missing Information)
# Formula: Mean(diff(Energy)^2) / Var(Energy)
bfmi_vector <- function(energy) {
  # Requires at least 3 samples to calculate diff and variance reliably
  if (length(energy) < 3) return(NA_real_)
  # Numerator: mean squared difference between consecutive energy levels
  num <- mean(diff(energy)^2)
  # Denominator: variance of energy levels
  den <- var(energy)
  if (is.na(den) || den == 0) return(NA_real_)
  num / den
}

bfmi_per_chain <- function(fit) {
  # Extract NUTS parameters (including energy__)
  nuts <- get_sampler_params(fit, inc_warmup = FALSE)
  # Apply bfmi_vector to the energy column of each chain's dataframe
  vapply(nuts, \(df) bfmi_vector(df[, "energy__"]), FUN.VALUE = numeric(1))
}

# -------------------------------------------------------------------------
# Main inspection routine
# -------------------------------------------------------------------------
inspect_fit <- function(path) {

  cat(">> Inspecting:", path, "\n")
  obj <- safe_read(path)

  # --- load errors -------------------------------------------------------
  if (inherits(obj, "error")) {
    cat("!! RDS read error:", obj$message, "\n"); return(invisible(NULL))
  }
  if (!inherits(obj, "stanfit")) {
    cat("!! Not a 'stanfit' object; classes:",
        paste(class(obj), collapse = ", "), "\n"); return(invisible(NULL))
  }
  # Check if the fit object is empty (e.g., compilation failed or 0 iterations)
  if (length(obj@sim) == 0 || is.null(obj@sim$iter) || obj@sim$iter == 0) {
    cat("!! stanfit contains zero iterations (empty fit)\n")
    return(invisible(NULL))
  }

  # --- basic meta --------------------------------------------------------
  # Use posterior package to extract draws and info
  draws   <- posterior::as_draws_array(obj)
  ndraw   <- posterior::ndraws(draws)
  nchains <- posterior::nchains(draws)
  n_iter  <- posterior::niterations(draws)

  cat("Model name            :", obj@model_name, "\n")
  cat("Chains × draws (post) :", nchains, "×", n_iter, "=", ndraw, "\n")

  # --- sampler diagnostics (NUTS) ----------------------------------------
  sp      <- get_sampler_params(obj, inc_warmup = FALSE)
  
  # 1. Divergences
  diverg  <- sum(vapply(sp, \(x) sum(x[, "divergent__"]), 0))
  
  # 2. Treedepth saturation
  # Get max_treedepth setting from the arguments used for the fit
  max_td_setting <- obj@stan_args[[1]]$control$max_treedepth %||% 10 # Default 10 if not found
  treedp  <- sum(vapply(sp, \(x) sum(x[, "treedepth__"] >= max_td_setting), 0))

  # 3. E-BFMI
  bfmi <- bfmi_per_chain(obj)
  
  cat("Divergent transitions :", diverg, paste0("(", round(100*diverg/ndraw, 2), "%)"), "\n")
  cat("Treedepth saturations :", treedp,   " (max =", max_td_setting, ")\n")
  cat("Energy BFMI per chain :", paste(round(bfmi, 3), collapse = ", "), "\n")
  if (any(bfmi < 0.2)) {
      cat("!! WARNING: Low E-BFMI (<0.2) detected (indicates poor exploration)\n")
  }

  # --- posterior summary (Convergence) -----------------------------------
  # Calculate Rhat and ESS using the posterior package (more robust methods)
  
  # Use summarize_draws for a comprehensive summary
  # Suppress warnings for parameters with 0 variance (e.g., fixed parameters)
  diags <- suppressWarnings(posterior::summarise_draws(draws, default_convergence_measures()))

  if (nrow(diags) == 0) {
      cat("!! Could not calculate convergence diagnostics.\n")
      return(invisible(list(fit = obj, draws = draws, sampler_params = sp)))
  }
  
  Rhat_max <- max(diags$rhat, na.rm = TRUE)
  cat("Worst R‑hat           :", round(Rhat_max, 3), "\n")

  min_bulk <- round(min(diags$ess_bulk, na.rm = TRUE))
  min_tail <- round(min(diags$ess_tail, na.rm = TRUE))
  
  cat("Min bulk ESS          :", min_bulk, "\n")
  cat("Min tail ESS          :", min_tail, "\n\n")

  # list top offenders by R‑hat (Rhat > 1.01 is a common threshold)
  bad_params <- filter(diags, rhat > 1.01)
  
  if (nrow(bad_params) > 0) {
    cat("Parameters with Rhat > 1.01:\n")
    out <- bad_params |>
      select(param = variable, Rhat = rhat, bulk_ESS = ess_bulk, tail_ESS = ess_tail) |>
      arrange(desc(Rhat))
    # Print without tibble truncation
    print(out, n = nrow(out))
  } else {
    cat("All R‑hat ≤ 1.01\n")
  }

  # Return the objects invisibly for further interactive analysis
  invisible(list(fit = obj, draws = draws, diagnostics = diags, sampler_params = sp))
}

# -------------------------------------------------------------------------
# Command‑line interface
# -------------------------------------------------------------------------
# This block executes if the script is run from the command line (Rscript)
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    cat("USAGE: Rscript inspect_fit.R path/to/fit_file.rds\n")
    # Exit with error status
    quit(status = 1)
  }
  # Check if file exists before attempting inspection
  if (!file.exists(args[1])) {
      cat("ERROR: File not found at", args[1], "\n")
      quit(status = 1)
  }
  invisible(inspect_fit(args[1]))
}

# -------------------------------------------------------------------------
# Interactive demo (only when sourced interactively)
# -------------------------------------------------------------------------
if (interactive()) {
  # Try to find an example file in the 'fits' directory for demo purposes
  if (dir.exists("fits")) {
      demo_fp <- list.files("fits", "^fit_.*\\.rds$", full.names = TRUE)[1]
      if (!is.na(demo_fp)) {
        cat("\n--- Running demo on", demo_fp, "---\n")
        inspect_fit(demo_fp)
      }
  }
}