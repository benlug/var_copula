#!/usr/bin/env Rscript
###########################################################################
# visualize_results.R  –  Differentiated by skewness + T × rho
#   • outcome metrics (all reps)
#   • outcome metrics (split by clean/problem)
#   • replication‑issue bar chart
#   all stored beneath  results/plots/<skew_level>/<T_rho>/...
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(rstan)
})

# ── Paths ----------------------------------------------------------------
BASE <- this.path::this.dir()
DIR <- list(
  data  = file.path(BASE, "data"),
  fits  = file.path(BASE, "fits"),
  res   = file.path(BASE, "results"),
  plots = file.path(BASE, "results", "plots")
)
dir.create(DIR$plots, FALSE, TRUE)

files <- list(
  cond   = file.path(DIR$res, "summary_conditions.csv"),
  rep    = file.path(DIR$res, "summary_replications.csv"),
  design = file.path(DIR$data, "sim_conditions_singlelevel.rds")
)
stopifnot(
  file.exists(files$cond), file.exists(files$rep),
  file.exists(files$design)
)

# ── Data -----------------------------------------------------------------
design <- readRDS(files$design) |>
  select(condition_id, skew_level, direction, T, rho, VARset)

cond <- read_csv(files$cond, show_col_types = FALSE) |>
  left_join(design, by = "condition_id")

rep <- read_csv(files$rep, show_col_types = FALSE) |>
  filter(!is.na(param)) |>
  left_join(design, by = "condition_id")

# ── Robust max‑Rhat flag -------------------------------------------------
max_rhat_bad <- function(path, thr = 1.01) {
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || !inherits(obj, "stanfit")) {
    return(NA)
  }
  it <- tryCatch(obj@sim$iter, error = function(e) NA_integer_)
  if (length(it) != 1 || is.na(it) || it == 0) {
    return(NA)
  }
  rhat_max <- suppressWarnings(
    max(summary(obj)$summary[, "Rhat"], na.rm = TRUE)
  )
  if (is.infinite(rhat_max)) {
    return(NA)
  }
  rhat_max > thr
}

# ── Append replication status -------------------------------------------
fit_files <- list.files(DIR$fits, "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$", TRUE)
rhat_tbl <- tibble(
  fit_path = fit_files,
  bad_rhat = map_lgl(fit_files, max_rhat_bad)
) |>
  mutate(
    model = str_match(basename(fit_path), "fit_(SG|NG)_")[, 2],
    condition_id = as.integer(str_match(fit_path, "cond(\\d+)")[, 2]),
    rep_id = as.integer(str_match(fit_path, "rep(\\d+)")[, 2])
  ) |>
  select(model, condition_id, rep_id, bad_rhat)

rep <- rep |>
  left_join(rhat_tbl, by = c("model", "condition_id", "rep_id")) |>
  mutate(
    bad_rhat  = coalesce(bad_rhat, FALSE),
    divergent = n_div > 0,
    status    = ifelse(bad_rhat | divergent, "problem", "clean")
  )

# ── Faceting helpers -----------------------------------------------------
facet_basic <- facet_grid(direction ~ VARset)
facet_status <- facet_grid(status + direction ~ VARset,
  labeller = labeller(
    status = c(clean = "clean", problem = "problem")
  )
)

plot_metric <- function(df, metric, ylab, out_file, facet_fun) {
  g <- ggplot(df, aes(param, .data[[metric]], fill = model)) +
    geom_col(position = position_dodge(width = .7), width = .6) +
    coord_flip() +
    facet_fun +
    theme_bw(base_size = 8) +
    labs(x = "", y = ylab, title = metric)
  ggsave(out_file, g, width = 11, height = 7)
}

# ── Loop structure: skew_level → T × rho --------------------------------
for (sk in unique(design$skew_level)) {
  sk_dir <- file.path(DIR$plots, sk)
  dir.create(sk_dir, FALSE, TRUE)

  cond_sk <- cond |> filter(skew_level == sk)
  rep_sk <- rep |> filter(skew_level == sk)

  combos <- unique(cond_sk[c("T", "rho")])

  for (i in seq_len(nrow(combos))) {
    Tval <- combos$T[i]
    rho_val <- combos$rho[i]
    tag <- paste0("T", Tval, "_rho", rho_val)
    out_d <- file.path(sk_dir, tag)
    dir.create(out_d, FALSE, TRUE)

    sel_c <- cond_sk |> filter(T == Tval, rho == rho_val)
    sel_r <- rep_sk |> filter(T == Tval, rho == rho_val)

    # -------- all reps --------------------------------------------------
    plot_metric(
      sel_c, "mean_rel_bias", "relative bias",
      file.path(out_d, "bias_rel.pdf"), facet_basic
    )
    plot_metric(
      sel_c, "coverage_95", "coverage",
      file.path(out_d, "coverage.pdf"), facet_basic
    )
    plot_metric(
      sel_c, "mean_post_sd", "posterior SD",
      file.path(out_d, "post_sd.pdf"), facet_basic
    )
    plot_metric(
      sel_c, "sd_bias", "SD‑Bias",
      file.path(out_d, "sd_bias.pdf"), facet_basic
    )
    plot_metric(
      sel_c, "mean_n_div", "# divergences",
      file.path(out_d, "divergences.pdf"), facet_basic
    )

    # -------- split by status ------------------------------------------
    rep_agg <- sel_r |>
      group_by(condition_id, model, param, status, direction, VARset) |>
      summarise(
        mean_rel_bias = mean(rel_bias),
        coverage_95 = mean(cover95),
        mean_post_sd = mean(post_sd),
        emp_sd = sd(post_mean),
        sd_bias = mean_post_sd - emp_sd,
        mean_n_div = mean(n_div),
        .groups = "drop"
      )

    plot_metric(
      rep_agg, "mean_rel_bias", "relative bias",
      file.path(out_d, "bias_rel_split.pdf"), facet_status
    )
    plot_metric(
      rep_agg, "coverage_95", "coverage",
      file.path(out_d, "coverage_split.pdf"), facet_status
    )
    plot_metric(
      rep_agg, "mean_post_sd", "posterior SD",
      file.path(out_d, "post_sd_split.pdf"), facet_status
    )
    plot_metric(
      rep_agg, "sd_bias", "SD‑Bias",
      file.path(out_d, "sd_bias_split.pdf"), facet_status
    )
    plot_metric(
      rep_agg, "mean_n_div", "# divergences",
      file.path(out_d, "divergences_split.pdf"), facet_status
    )
  }
}

message("✓ Plots saved under ", DIR$plots)
