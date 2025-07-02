#!/usr/bin/env Rscript
###########################################################################
# inspect_fit.R  –  Inspect a single Stan .rds fit file
#
# Usage:
#   Rscript inspect_fit.R path/to/fit_file.rds
# or interactively
#   source("inspect_fit.R"); inspect_fit("fits/fit_XXX.rds")
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tibble)
  library(posterior)
})

safe_read <- function(path) tryCatch(readRDS(path), error = function(e) e)

# ------------ main inspection function ----------------------------------
inspect_fit <- function(path) {
  cat(">> Inspecting:", path, "\n")
  obj <- safe_read(path)

  ## --- error at load ---------------------------------------------------
  if (inherits(obj, "error")) {
    cat("!! Could not read file (RDS error):", obj$message, "\n")
    return(invisible(NULL))
  }

  ## --- not a stanfit ---------------------------------------------------
  if (!inherits(obj, "stanfit")) {
    cat(
      "!! Object is not of class 'stanfit' (class:",
      paste(class(obj), collapse = ", "), ")\n"
    )
    return(invisible(NULL))
  }

  ## --- empty draw set --------------------------------------------------
  if (obj@sim$iter == 0) {
    cat("!! stanfit contains zero iterations (fit skipped or failed)\n")
    return(invisible(NULL))
  }

  ## --- basic meta ------------------------------------------------------
  draws <- posterior::as_draws_array(obj)
  ndraw <- posterior::ndraws(draws)
  nchains <- posterior::nchains(draws)
  cat("Model name     :", obj@model_name, "\n")
  cat(
    "Chains × draws :", nchains, "×", ndraw / nchains,
    "(total", ndraw, ")\n"
  )

  ## --- sampler params --------------------------------------------------
  sp <- rstan::get_sampler_params(obj, FALSE)
  diverg <- sum(vapply(sp, \(x) sum(x[, "divergent__"]), 0))
  treedp <- sum(vapply(sp, \(x) sum(x[, "treedepth__"] >= obj@stan_args[[1]]$control$max_treedepth), 0))

  cat("Divergent transitions :", diverg, "\n")
  cat("Treedepth saturations :", treedp, "\n")

  ## energy BFMI
  energies <- unlist(lapply(sp, \(x) x[, "energy__"]))
  bfmi <- if (length(energies) > 1) {
    stats::var(energies) / mean(diff(energies)^2)
  } else {
    NA_real_
  }
  cat("Energy BFMI            :", round(bfmi, 3), "\n")

  ## --- summary stats ---------------------------------------------------
  summ <- summary(obj)$summary
  worst_rhat <- head(sort(summ[, "Rhat"], decreasing = TRUE), 10)
  worst_rhat <- worst_rhat[!is.na(worst_rhat)]

  bulk <- posterior::ess_bulk(draws)
  tail <- posterior::ess_tail(draws)

  cat(
    "Max R-hat              :", max(summ[, "Rhat"], na.rm = TRUE),
    "\n"
  )
  cat("Min bulk ESS           :", min(bulk), "\n")
  cat("Min tail ESS           :", min(tail), "\n\n")

  if (length(worst_rhat)) {
    cat("Top parameters by R-hat (>1.01):\n")
    idx <- which(summ[, "Rhat"] %in% worst_rhat)
    out <- tibble(
      param = rownames(summ)[idx],
      Rhat = summ[idx, "Rhat"],
      bulk_ESS = bulk[idx],
      tail_ESS = tail[idx]
    ) |>
      arrange(desc(Rhat))
    print(out, n = nrow(out))
  } else {
    cat("All R-hat ≤ 1.01\n")
  }

  invisible(list(
    fit       = obj,
    sampler   = sp,
    summary   = summ
  ))
}

# ------------ command‑line interface ------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("USAGE: Rscript inspect_fit.R path/to/fit_file.rds\n")
    quit(status = 1)
  }
  invisible(inspect_fit(args[1]))
}


BASE_DIR <- this.path::this.dir()
setwd(BASE_DIR)
getwd()
res <- inspect_fit(paste0(BASE_DIR, "/fits/fit_SG_cond001_rep001.rds"))
res$summary
