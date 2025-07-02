###########################################################################
# diagnose_fits.R  –  Summarise Stan diagnostics for every .rds fit
#
# * Collects: divergences, R-hat, bulk/tail ESS
# * Writes results/diagnostics_summary.csv  (one row per fit)
# * Run after fit_models.R has created files in fits/
###########################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(stringr)
  library(readr)
})

# --------------- user-editable directories ------------------------------
BASE_DIR <- this.path::this.dir() # script location
FITS_DIR <- file.path(BASE_DIR, "fits")
RES_DIR <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, FALSE, TRUE)

# --------------- helper: count divergent transitions --------------------
count_divergent <- function(fit) {
  sum(vapply(
    rstan::get_sampler_params(fit, inc_warmup = FALSE),
    function(x) sum(x[, "divergent__"]), numeric(1)
  ))
}

# --------------- iterate over fit files ---------------------------------
fit_files <- list.files(FITS_DIR,
  "^fit_(SNG|NG)_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)

if (!length(fit_files)) {
  stop("No Stan fit files found in ", FITS_DIR)
}

diagnostics <- purrr::map_dfr(fit_files, function(ff) {
  meta <- str_match(
    basename(ff),
    "fit_(SNG|NG)_cond(\\d+)_rep(\\d+)\\.rds"
  )
  model_code <- meta[2]
  cond_id <- as.integer(meta[3])
  rep_id <- as.integer(meta[4])

  fit <- readRDS(ff)

  sum_df <- as.data.frame(summary(fit)$summary)
  rhat <- sum_df$Rhat
  ess <- sum_df$n_eff # fallback for pre‑2.26 rstan
  if (!is.null(sum_df$ess_bulk__)) {
    ess_bulk <- sum_df$ess_bulk__
    ess_tail <- sum_df$ess_tail__
  } else {
    ess_bulk <- ess_tail <- ess # older Stan versions
  }

  tibble(
    model        = model_code,
    condition_id = cond_id,
    rep_id       = rep_id,
    n_divergent  = count_divergent(fit),
    max_rhat     = max(rhat, na.rm = TRUE),
    pct_rhat_gt  = mean(rhat > 1.01, na.rm = TRUE),
    min_bulk_ess = min(ess_bulk, na.rm = TRUE),
    med_bulk_ess = median(ess_bulk, na.rm = TRUE),
    min_tail_ess = min(ess_tail, na.rm = TRUE),
    med_tail_ess = median(ess_tail, na.rm = TRUE)
  )
})

out_file <- file.path(RES_DIR, "diagnostics_summary.csv")
write_csv(diagnostics, out_file)
message("✓ Diagnostics written: ", out_file)
