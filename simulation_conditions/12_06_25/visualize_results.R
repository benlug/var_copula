#!/usr/bin/env Rscript
###########################################################################
# visualize_results.R  –  Updated July‑2025
#   – works with new σ‑parameters
#   – one PDF page per figure
#   – tables saved alongside plots
#   – global (across‑T) trend plots + matching box‑plots
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

## ── paths ---------------------------------------------------------------
BASE <- this.path::this.dir()
DIR <- list(
  data = file.path(BASE, "data"),
  fits = file.path(BASE, "fits"),
  res = file.path(BASE, "results"),
  plots = file.path(BASE, "results", "plots"),
  plots_g = file.path(BASE, "results", "plots_global")
)
dir.create(DIR$plots, FALSE, TRUE)
dir.create(DIR$plots_g, FALSE, TRUE)

files <- list(
  cond   = file.path(DIR$res, "summary_conditions.csv"),
  rep    = file.path(DIR$res, "summary_replications.csv"),
  design = file.path(DIR$data, "sim_conditions_singlelevel.rds")
)
stopifnot(
  file.exists(files$cond), file.exists(files$rep),
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

## ── parameter ordering --------------------------------------------------
param_levels <- c(
  "omega[1]", "omega[2]", "alpha[1]", "alpha[2]",
  "sigma[1]", "sigma[2]",
  "mu[1]", "mu[2]",
  "phi11", "phi12", "phi21", "phi22",
  "rho"
)
rep <- mutate(rep, param = factor(param, levels = param_levels))
cond <- mutate(cond, param = factor(param, levels = param_levels))

## ── divergence & R‑hat status ------------------------------------------
max_rhat_bad <- function(path, thr = 1.01) {
  if (!file.exists(path)) {
    return(NA)
  }
  o <- tryCatch(readRDS(path), error = \(e) NULL)
  if (!inherits(o, "stanfit")) {
    return(NA)
  }
  it <- tryCatch(o@sim$iter, error = \(e) NA_integer_)
  if (is.na(it) || it == 0) {
    return(NA)
  }
  rmax <- suppressWarnings(max(summary(o)$summary[, "Rhat"], na.rm = TRUE))
  if (is.infinite(rmax)) {
    return(NA)
  }
  rmax > thr
}
ff <- list.files(DIR$fits, "^fit_(SG|NG)_cond\\d+_rep\\d+\\.rds$", TRUE)
rhat_tbl <- tibble(
  fit = ff,
  bad_rhat = map_lgl(ff, max_rhat_bad)
) |>
  mutate(
    model = str_match(basename(fit), "fit_(SG|NG)_")[, 2],
    condition_id = as.integer(str_match(fit, "cond(\\d+)")[, 2]),
    rep_id = as.integer(str_match(fit, "rep(\\d+)")[, 2])
  ) |>
  select(-fit)

rep <- rep |>
  left_join(rhat_tbl, by = c("model", "condition_id", "rep_id")) |>
  mutate(
    bad_rhat = coalesce(bad_rhat, FALSE),
    divergent = n_div > 0,
    status = ifelse(bad_rhat | divergent, "problem", "clean")
  )

## ── plotting helpers ----------------------------------------------------
facet_basic <- facet_grid(direction ~ VARset,
  labeller = labeller(VARset = label_both)
)
facet_status <- facet_grid(status + direction ~ VARset,
  labeller = labeller(
    VARset = label_both,
    status = c(
      clean = "clean",
      problem = "problem"
    )
  )
)

tab_grob <- function(df) {
  gridExtra::tableGrob(df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = grid::unit(c(2, 2), "mm")
    )
  )
}

write_metric_page <- function(df, metric, ylab, file_pdf,
                              facet_fun, with_table = FALSE) {
  g <- ggplot(df, aes(param, .data[[metric]], fill = model)) +
    geom_col(position = position_dodge(0.65), width = .55) +
    coord_flip() +
    facet_fun +
    theme_bw(8) +
    labs(x = "", y = ylab, title = metric)

  pdf(file_pdf, 11, 7)
  if (with_table) {
    gridExtra::grid.arrange(
      tab_grob(df |> select(
        condition_id, model, param,
        !!metric := .data[[metric]]
      ))
    )
  }
  print(g)
  dev.off()
}

write_metric_split <- function(rep_df, metric, ylab, file_pdf) {
  rep_df <- mutate(rep_df,
    .val = case_when(
      metric == "rel_bias" ~ rel_bias,
      metric == "coverage" ~ cover95,
      metric == "post_sd" ~ post_sd,
      metric == "sd_bias" ~ post_sd, # temp, adjust later
      metric == "n_div" ~ n_div,
      TRUE ~ NA_real_
    )
  )

  if (metric == "sd_bias") {
    rep_df <- rep_df |>
      group_by(condition_id, model, param) |>
      mutate(
        emp_sd = sd(post_mean),
        .val = post_sd - emp_sd
      ) |>
      ungroup()
  }

  tbl_cnt <- rep_df |> count(status, direction, VARset, name = "n_reps")

  agg <- rep_df |>
    group_by(condition_id, model, param, status, direction, VARset) |>
    summarise(mean_val = mean(.val, na.rm = TRUE), .groups = "drop")

  pdf(file_pdf, 11, 7)
  gridExtra::grid.arrange(tab_grob(tbl_cnt))
  print(
    ggplot(agg, aes(param, mean_val, fill = model)) +
      geom_col(position = "dodge", width = .55) +
      coord_flip() +
      facet_status +
      theme_bw(8) +
      labs(
        x = "", y = ylab,
        title = paste0(metric, " (clean vs problem)")
      )
  )
  dev.off()
}

write_div_table <- function(rep_df, file_pdf) {
  tbl <- rep_df |>
    group_by(condition_id, model, direction, VARset) |>
    summarise(n_div = sum(n_div), n_reps = n(), .groups = "drop")
  pdf(file_pdf, 11, 7)
  gridExtra::grid.arrange(tab_grob(tbl))
  dev.off()
}

global_plots <- function(df, metric, ylab, sk_tag) {
  for (rho_val in sort(unique(df$rho))) {
    df_r <- df |> filter(rho == rho_val)
    tag <- paste0(sk_tag, "_rho", rho_val, "_", metric)
    ## line
    g1 <- ggplot(df_r, aes(factor(T), .data[[metric]],
      colour = model, group = model
    )) +
      geom_line() +
      geom_point() +
      facet_grid(param ~ direction + VARset,
        scales = "free_y",
        labeller = labeller(VARset = label_both)
      ) +
      theme_bw(7) +
      labs(
        x = "T", y = ylab,
        title = paste0(metric, "  |  ρ=", rho_val)
      )
    ## box
    g2 <- ggplot(df_r, aes(factor(T), .data[[metric]], fill = model)) +
      geom_boxplot(outlier.size = .4, position = "dodge") +
      facet_grid(param ~ direction + VARset,
        scales = "free_y",
        labeller = labeller(VARset = label_both)
      ) +
      theme_bw(7) +
      labs(
        x = "T", y = ylab,
        title = paste0("box‑plot ", metric, "  |  ρ=", rho_val)
      )

    pdf(file.path(DIR$plots_g, paste0(tag, ".pdf")), 12, 8)
    print(g1)
    print(g2)
    dev.off()
  }
}

## ── glossary ------------------------------------------------------------
writeLines(c(
  "Metric definitions:",
  "  relative bias  = (posterior mean − truth) / |truth|",
  "  coverage_95    = proportion of reps where 95% CI covers truth",
  "  posterior SD   = mean posterior standard deviation across reps",
  "  SD‑Bias        = posterior SD − empirical SD(post. means)",
  "  n_div          = divergent transitions (after warm‑up)",
  "  clean          = no divergence & Rhat ≤ 1.01"
), file.path(DIR$res, "metric_glossary.txt"))

## ── main loop -----------------------------------------------------------
metrics_c <- list(
  mean_rel_bias = list(ylab = "relative bias", m = "rel_bias"),
  coverage_95   = list(ylab = "coverage", m = "coverage"),
  mean_post_sd  = list(ylab = "posterior SD", m = "post_sd"),
  sd_bias       = list(ylab = "SD‑Bias", m = "sd_bias"),
  mean_n_div    = list(ylab = "# divergences", m = "n_div")
)

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
    outd <- file.path(sk_dir, tag)
    dir.create(outd, FALSE, TRUE)

    sel_c <- cond_sk |> filter(T == Tval, rho == rho_val)
    sel_r <- rep_sk |> filter(T == Tval, rho == rho_val)

    ## unsplit pages + tables
    lapply(names(metrics_c), function(mc) {
      write_metric_page(
        sel_c, mc, metrics_c[[mc]]$ylab,
        file.path(outd, paste0(mc, ".pdf")),
        facet_basic,
        with_table = TRUE
      )
    })

    ## split pages
    lapply(names(metrics_c), function(mc) {
      write_metric_split(
        sel_r, metrics_c[[mc]]$m, metrics_c[[mc]]$ylab,
        file.path(outd, paste0(mc, "_split.pdf"))
      )
    })

    ## divergence count table
    write_div_table(sel_r, file.path(outd, "divergence_table.pdf"))
  }

  ## global trend plots
  for (mc in names(metrics_c)) {
    global_plots(cond_sk, mc, metrics_c[[mc]]$ylab, sk)
  }
}

message(
  "✓ visualisation complete – see ",
  DIR$plots, " and ", DIR$plots_g
)
