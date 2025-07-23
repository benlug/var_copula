## quick_fit_check.R ------------------------------------------------------
##
## Safely open a fit_*.rds file and print:
##   • status  (ok / empty / corrupt)
##   • Stan error message if available
## -----------------------------------------------------------------------

inspect_fit_file_once <- function(path) {
  cat("Opening:", path, "\n")

  if (!file.exists(path)) {
    cat("  -> FILE MISSING\n")
    return(invisible(NULL))
  }

  obj <- tryCatch(readRDS(path), error = function(e) e)

  if (inherits(obj, "error")) {
    cat("  -> CORRUPT .rds (readRDS error)\n")
    cat("  ", conditionMessage(obj), "\n")
    return(invisible(NULL))
  }

  if (!inherits(obj, "stanfit")) {
    cat("  -> CORRUPT .rds (not a stanfit object)\n")
    cat("  class(obj) = ", paste(class(obj), collapse = ", "), "\n")
    return(invisible(NULL))
  }

  it <- tryCatch(obj@sim$iter, error = function(e) NA_integer_)
  if (length(it) != 1 || is.na(it) || it == 0) {
    cat("  -> EMPTY stanfit (0 iterations)\n")
    # Stan’s low‑level error message (cmdstan/rstan >=2.32) is stored here:
    msg <- tryCatch(obj@stan_args[[1]]$msg, error = function(e) NULL)
    if (!is.null(msg)) cat("  Stan message:\n  ", msg, "\n", sep = "")
    return(invisible(NULL))
  }

  n_div <- sum(vapply(
    rstan::get_sampler_params(obj, FALSE),
    function(x) sum(x[, "divergent__"]), 0
  ))
  cat("  -> OK  (", it, " iterations ; ", n_div, " divergences)\n", sep = "")
  invisible(NULL)
}

BASE_DIR <- this.path::this.dir()
path <- file.path(BASE_DIR, "../fits/fit_SG_cond001_rep008.rds")
inspect_fit_file_once(path)

fit@stan_args[[k]]$seed


# =========================================================================
# safe_extract_seed_init()      – one‑stop wrapper
# -------------------------------------------------------------------------
# • Returns a list( ok = TRUE/FALSE,
#                  seed = <int or NA>,
#                  init = <list | NULL>,
#                  msg  = explanatory character )
# • Handles all corner cases found in your fits/ folder.
# =========================================================================

safe_extract_seed_init <- function(path, chain = 1) {
  if (!file.exists(path)) {
    return(list(
      ok = FALSE, seed = NA, init = NULL,
      msg = "file does not exist"
    ))
  }

  obj <- tryCatch(readRDS(path), error = function(e) e)

  ## ---- unreadable or not an R object -----------------------------------
  if (inherits(obj, "error")) {
    return(list(
      ok = FALSE, seed = NA, init = NULL,
      msg = paste("readRDS error:", conditionMessage(obj))
    ))
  }

  ## ---- some wrappers save a *character* in place of a stanfit ----------
  if (!inherits(obj, "stanfit")) {
    return(list(
      ok = FALSE, seed = NA, init = NULL,
      msg = paste(
        "object class is", paste(class(obj), collapse = ", "),
        "– not a stanfit"
      )
    ))
  }

  ## ---- deal with empty stanfit: stan_args may be length‑0 --------------
  if (length(obj@stan_args) < chain) {
    return(list(
      ok = FALSE, seed = NA, init = NULL,
      msg = "stanfit present but stan_args slot is empty (0‑iter run)"
    ))
  }

  sarg <- obj@stan_args[[chain]]
  seed <- sarg$seed %||% NA_integer_

  init_obj <- sarg$init
  init_list <- NULL

  if (is.function(init_obj)) {
    if (!is.na(seed)) {
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    }
    init_list <- tryCatch(init_obj(), error = function(e) NULL)
  } else if (!is.null(init_obj)) {
    init_list <- init_obj
  }

  if (is.null(init_list)) {
    return(list(
      ok = FALSE, seed = seed, init = NULL,
      msg = "could not materialise init (function errored?)"
    ))
  }

  list(
    ok = TRUE, seed = seed, init = init_list,
    msg = "seed & init extracted"
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b # tiny infix helper


BASE_DIR <- this.path::this.dir() # folder that contains the script
path <- file.path(
  BASE_DIR,
  "../fits/fit_SG_cond001_rep008.rds"
)

info <- safe_extract_seed_init(path, chain = 1)

if (info$ok) {
  print(info$seed)
  str(info$init)
} else {
  cat("Could not extract seed/init:\n", info$msg, "\n")
}
