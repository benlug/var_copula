###########################################################################
# investigate_test_fits.R
#
# Inspect Stan fits produced by fit_test_data.R.
# Creates summaries of convergence diagnostics and parameter estimates
# and saves simple traceplots for each parameter.
###########################################################################


# --- Libraries ---
library(rstan)
library(dplyr)
library(ggplot2)

# Determine script directory and base directory for results
SCRIPT_DIR <- tryCatch(this.path::this.dir(), error = function(e) getwd())

# --- Configuration ---
fits_dir <- file.path(SCRIPT_DIR, "test_fits")
plots_dir <- file.path(SCRIPT_DIR, "fit_diagnostics")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

fit_files <- list.files(fits_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(fit_files) == 0) {
  stop("No fit files found in ", fits_dir)
}

# container for summaries
summaries <- list()

# --- Helper for plotting ---
make_traceplot <- function(fit, param, file_stub) {
  p <- rstan::traceplot(fit, pars = param, inc_warmup = FALSE) +
    ggtitle(sprintf("Traceplot for %s", param))
  ggsave(file.path(plots_dir, paste0(file_stub, "_", param, "_trace.png")),
    plot = p, width = 6, height = 4
  )
}

# --- Iterate over fits ---
for (ff in fit_files) {
  fit <- readRDS(ff)
  summ <- summary(fit)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "parameter") %>%
    dplyr::select(parameter, mean, sd, n_eff = n_eff, Rhat)
  summaries[[basename(ff)]] <- summ

  # make traceplots for all parameters except lp__
  params <- setdiff(summ$parameter, "lp__")
  for (param in params) {
    file_stub <- tools::file_path_sans_ext(basename(ff))
    make_traceplot(fit, param, file_stub)
  }
}

# combine and save summary table
summary_df <- dplyr::bind_rows(summaries, .id = "fit_file")
summary_out <- file.path(plots_dir, "fit_summaries.csv")
readr::write_csv(summary_df, summary_out)
cat("Saved summary diagnostics to", summary_out, "\n")
