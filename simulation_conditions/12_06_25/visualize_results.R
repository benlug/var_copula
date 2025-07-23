#!/usr/bin/env Rscript
###########################################################################
# visualize_results.R  –  robust version (no patchwork, no sd_bias bug)
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(gridExtra)
  library(rstan)
})

## ── project paths ──────────────────────────────────────────────────────
BASE <- this.path::this.dir()
DIR <- list(
  data = file.path(BASE, "data"),
  fits = file.path(BASE, "fits"),
  res = file.path(BASE, "results"),
  plots = file.path(BASE, "results", "plots"),
  plots_global = file.path(BASE, "results", "plots_global")
)
dir.create(DIR$plots, FALSE, TRUE)
dir.create(DIR$plots_global, FALSE, TRUE)

files <- list(
  cond   = file.path(DIR$res, "summary_conditions.csv"),
  rep    = file.path(DIR$res, "summary_replications.csv"),
  design = file.path(DIR$data, "sim_conditions_singlelevel.rds")
)
stopifnot(
  file.exists(files$cond),
  file.exists(files$rep),
  file.exists(files$design)
)

## ── data ----------------------------------------------------------------
design <- readRDS(files$design) |>
  select(condition_id, skew_level, direction, T, rho, VARset)

cond <- read_csv(files$cond, show_col_types = FALSE) |>
  left_join(design, by = "condition_id")

rep <- read_csv(files$rep, show_col_types = FALSE) |>
  filter(!is.na(param)) |>
  left_join(design, by = "condition_id")

## ── append status flags (divergences + bad Rhat) ------------------------
max_rhat_bad <- function(path, thr = 1.01) {
  if (!file.exists(path)) {
    return(NA)
  }
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || !inherits(obj, "stanfit")) {
    return(NA)
  }
  it <- tryCatch(obj@sim$iter, error = function(e) NA_integer_)
  if (length(it) != 1 || is.na(it) || it == 0) {
    return(NA)
  }
  rh <- suppressWarnings(max(summary(obj)$summary[, "Rhat"], na.rm = TRUE))
  if (is.infinite(rh)) {
    return(NA)
  }
  rh > thr
}

fit_paths <- list.files(DIR$fits,
  "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)

rhat_tbl <- tibble(
  fit_path = fit_paths,
  bad_rhat = map_lgl(fit_paths, max_rhat_bad)
) |>
  mutate(
    model        = str_match(basename(fit_path), "fit_(SG|NG)_")[, 2],
    condition_id = as.integer(str_match(fit_path, "cond(\\d+)")[, 2]),
    rep_id       = as.integer(str_match(fit_path, "rep(\\d+)")[, 2])
  ) |>
  select(-fit_path)

rep <- rep |>
  left_join(rhat_tbl, by = c("model", "condition_id", "rep_id")) |>
  mutate(
    bad_rhat  = coalesce(bad_rhat, FALSE),
    divergent = n_div > 0,
    status    = ifelse(bad_rhat | divergent, "problem", "clean")
  )

## ── plotting helpers ----------------------------------------------------
facet_basic <- facet_grid(direction ~ VARset)
facet_status <- facet_grid(status + direction ~ VARset,
  labeller = labeller(
    status = c(clean = "clean", problem = "problem")
  )
)

table_grob <- function(df) {
  gridExtra::tableGrob(df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = unit(c(2, 2), "mm")
    )
  )
}

write_metric_plot <- function(df, metric, ylab, out_file, facet_fun) {
  g <- ggplot(df, aes(param, .data[[metric]], fill = model)) +
    geom_col(position = position_dodge(0.7), width = .6) +
    coord_flip() +
    facet_fun +
    theme_bw(8) +
    labs(x = "", y = ylab, title = metric)
  ggsave(out_file, g, width = 11, height = 7)
}

## ---- split‑metric (table + plot) ---------------------------------------
write_metric_split <- function(rep_df, metric, ylab,
                               out_file, facet_fun) {
  ## --- compute metric at replication level -----------------------------
  rep_df <- rep_df |>
    mutate(.val = case_when(
      metric == "rel_bias" ~ rel_bias,
      metric == "cover95" ~ cover95,
      metric == "post_sd" ~ post_sd,
      metric == "n_div" ~ n_div,
      metric == "sd_bias" ~ NA_real_, # placeholder
      TRUE ~ NA_real_
    ))

  ## sd_bias needs posterior‑SD & empirical‑SD across reps: do it here
  if (metric == "sd_bias") {
    rep_df <- rep_df |>
      group_by(condition_id, model, param) |>
      mutate(
        emp_sd = sd(post_mean),
        sd_bias = post_sd - emp_sd
      ) |>
      ungroup() |>
      mutate(.val = sd_bias)
  }

  ## --- counts table -----------------------------------------------------
  cnt_tbl <- rep_df |>
    count(status, direction, VARset, name = "n_reps") |>
    arrange(VARset, direction, status)

  ## --- aggregation for bar‑plot ----------------------------------------
  rep_agg <- rep_df |>
    group_by(
      condition_id, model, param,
      status, direction, VARset
    ) |>
    summarise(mean_val = mean(.val, na.rm = TRUE), .groups = "drop")

  ## --- figure -----------------------------------------------------------
  g_tab <- table_grob(cnt_tbl)

  g_plot <- ggplot(
    rep_agg,
    aes(param, mean_val, fill = model)
  ) +
    geom_col(position = "dodge", width = .65) +
    coord_flip() +
    facet_fun +
    theme_bw(8) +
    labs(
      x = "", y = ylab,
      title = paste0(metric, "  (clean vs problem)")
    )

  pdf(out_file, 11, 7)
  gridExtra::grid.arrange(g_tab) # page 1
  print(g_plot) # page 2
  dev.off()
}

write_divergence_table <- function(rep_df, out_file) {
  tbl <- rep_df |>
    group_by(condition_id, model, direction, VARset) |>
    summarise(
      n_divergent = sum(n_div),
      n_reps = n(), .groups = "drop"
    )
  pdf(out_file, 11, 7)
  gridExtra::grid.arrange(table_grob(tbl))
  dev.off()
}

## ── metric glossary -----------------------------------------------------
writeLines(c(
  "Glossary of Simulation Metrics",
  "-------------------------------------------",
  "relative bias  = (posterior mean – truth) / |truth|",
  "coverage_95    = 1 if 95% CI encloses truth, else 0",
  "posterior SD   = mean posterior standard deviation",
  "SD‑Bias        = posterior SD – empirical SD across reps",
  "n_div          = # divergent transitions (after warm‑up)",
  "status         = 'clean' (no divergences & good Rhat) vs problem"
), file.path(DIR$res, "metric_glossary.txt"))

## ── per‑skew‑level loop --------------------------------------------------
for (sk in unique(design$skew_level)) {
  sk_dir <- file.path(DIR$plots, sk)
  dir.create(sk_dir, FALSE, TRUE)
  cond_sk <- cond |> filter(skew_level == sk)
  rep_sk <- rep |> filter(skew_level == sk)

  combos <- unique(cond_sk[c("T", "rho")])

  for (j in seq_len(nrow(combos))) {
    Tval <- combos$T[j]
    rho_val <- combos$rho[j]
    tag <- paste0("T", Tval, "_rho", rho_val)
    out_d <- file.path(sk_dir, tag)
    dir.create(out_d, FALSE, TRUE)

    sel_c <- cond_sk |> filter(T == Tval, rho == rho_val)
    sel_r <- rep_sk |> filter(T == Tval, rho == rho_val)

    ## --- simple bar‑plots (all reps) -----------------------------------
    write_metric_plot(
      sel_c, "mean_rel_bias", "relative bias",
      file.path(out_d, "bias_rel.pdf"), facet_basic
    )
    write_metric_plot(
      sel_c, "coverage_95", "coverage",
      file.path(out_d, "coverage.pdf"), facet_basic
    )
    write_metric_plot(
      sel_c, "mean_post_sd", "posterior SD",
      file.path(out_d, "post_sd.pdf"), facet_basic
    )
    write_metric_plot(
      sel_c, "sd_bias", "SD‑Bias",
      file.path(out_d, "sd_bias.pdf"), facet_basic
    )
    write_metric_plot(
      sel_c, "mean_n_div", "# divergences",
      file.path(out_d, "divergences.pdf"), facet_basic
    )

    ## --- split plots (table + bar) -------------------------------------
    write_metric_split(
      sel_r, "rel_bias", "relative bias",
      file.path(out_d, "bias_rel_split.pdf"), facet_status
    )
    write_metric_split(
      sel_r, "cover95", "coverage share",
      file.path(out_d, "coverage_split.pdf"), facet_status
    )
    write_metric_split(
      sel_r, "post_sd", "posterior SD",
      file.path(out_d, "post_sd_split.pdf"), facet_status
    )
    write_metric_split(
      sel_r, "sd_bias", "SD‑Bias",
      file.path(out_d, "sd_bias_split.pdf"), facet_status
    )
    write_metric_split(
      sel_r, "n_div", "# divergences",
      file.path(out_d, "divergences_split.pdf"), facet_status
    )

    ## --- stand‑alone divergence table ----------------------------------
    write_divergence_table(
      sel_r,
      file.path(out_d, "divergence_table.pdf")
    )
  }
}

## ── global (metric × T × ρ) plots ---------------------------------------
plot_global_metric <- function(df, metric, ylab, file_tag) {
  g <- ggplot(
    df,
    aes(interaction(T, rho, sep = "\nρ="),
      .data[[metric]],
      colour = model, group = model
    )
  ) +
    geom_point(position = position_dodge(0.3)) +
    geom_line(position = position_dodge(0.3)) +
    facet_grid(param ~ direction + VARset, scales = "free_y") +
    theme_bw(7) +
    labs(x = "T  × ρ", y = ylab, title = metric)
  ggsave(file.path(DIR$plots_global, paste0(file_tag, ".pdf")),
    g,
    width = 12, height = 8
  )
}

for (sk in unique(design$skew_level)) {
  df_sk <- cond |> filter(skew_level == sk)
  plot_global_metric(
    df_sk, "mean_rel_bias", "relative bias",
    paste0("global_bias_", sk)
  )
  plot_global_metric(
    df_sk, "coverage_95", "coverage",
    paste0("global_coverage_", sk)
  )
  plot_global_metric(
    df_sk, "mean_post_sd", "posterior SD",
    paste0("global_postsd_", sk)
  )
  plot_global_metric(
    df_sk, "sd_bias", "SD‑Bias",
    paste0("global_sdbias_", sk)
  )
  plot_global_metric(
    df_sk, "mean_n_div", "# divergences",
    paste0("global_ndiv_", sk)
  )
}

message(
  "✓ all plots & tables saved in ",
  DIR$plots, "  and  ", DIR$plots_global
)
