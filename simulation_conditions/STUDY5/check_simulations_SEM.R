# check_simulations_SEM.R — Visual checks directly on the simulated truth
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(patchwork)
})

# PIT for signed standardized exponential (scale=1)
exp_cdf_dir <- function(e, dir = +1L) {
  if (dir == +1L) { ifelse(e <= -1, 0, 1 - exp(-(e + 1))) }
  else            { ifelse(e >=  1, 1, exp(e - 1)) }
}
nscore <- function(u) stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))

theme_sm <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p_hist <- function(x, ttl, sub = "") ggplot(data.frame(x), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = .8) +
  labs(title = ttl, subtitle = sub, x = NULL, y = "density") + theme_sm

p_qq <- function(e, dir, ttl) {
  if (dir == -1L) e <- -e
  p <- (seq_along(e) - 0.5) / length(e)
  q <- -1 - log(1 - p)  # std Exp(+)
  df <- data.frame(theoretical = q, sample = sort(e))
  ggplot(df, aes(theoretical, sample)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = .6, size = 1.5) +
    labs(title = ttl, x = "Theoretical (std Exp+)", y = "Sample (mirrored if needed)") +
    theme_sm
}

p_pit <- function(u, ttl) ggplot(data.frame(u), aes(u)) +
  geom_histogram(bins = 20, color = "white") +
  geom_hline(yintercept = length(u) / 20, linetype = 2) +
  labs(title = ttl, x = "PIT(u)", y = "count") + theme_sm

p_cop <- function(u1, u2, rho, ttl) {
  z1 <- nscore(u1); z2 <- nscore(u2)
  ggplot(data.frame(z1, z2), aes(z1, z2)) +
    geom_point(alpha = .35, size = 1) +
    labs(title = ttl,
         subtitle = sprintf("empirical corr ≈ %.2f (true ρ=%.2f)", suppressWarnings(cor(z1, z2)), rho),
         x = "Φ^{-1}(u1)", y = "Φ^{-1}(u2)") +
    theme_sm
}

run_visual_checks_sem <- function(DATA_DIR = "data",
                                  OUT_DIR  = "checks_sem",
                                  reps_to_plot = 3) {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  paths <- list.files(DATA_DIR, "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(paths)) { message("No simulations in ", DATA_DIR); return(invisible()) }

  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  df   <- data.frame(path = paths, cid = as.integer(meta[,2]), rid = as.integer(meta[,3])) |>
          arrange(cid, rid) |>
          group_by(cid) |>
          slice_head(n = reps_to_plot) |>
          ungroup()

  for (k in seq_len(nrow(df))) {
    sim <- readRDS(df$path[k])
    cid <- df$cid[k]; rid <- df$rid[k]
    out_dir_cond <- file.path(OUT_DIR, sprintf("cond%03d", cid))
    dir.create(out_dir_cond, showWarnings = FALSE, recursive = TRUE)

    y   <- as.matrix(sim$data[, c("y1","y2")])
    dir_sign <- as.integer(sim$skew_signs)
    rho <- sim$rho
    Tt  <- nrow(y)

    if (identical(sim$sem_study, "A_indicator")) {
      e  <- sim$truth$meas_error
      u1 <- exp_cdf_dir(e[,1], dir_sign[1]); u2 <- exp_cdf_dir(e[,2], dir_sign[2])

      g <- ((p_hist(e[,1], "A: measurement residual ε1", sprintf("cond%03d rep%03d", cid, rid)) |
             p_hist(e[,2], "A: measurement residual ε2", sprintf("cond%03d rep%03d", cid, rid))) /
            (p_qq(e[,1], dir_sign[1], "QQ ε1 vs Exp±") |
             p_qq(e[,2], dir_sign[2], "QQ ε2 vs Exp±")) /
            (p_pit(u1, "PIT(u1)") | p_pit(u2, "PIT(u2)")) /
             p_cop(u1, u2, rho, "Gaussian-score scatter (measurement layer)")) +
        plot_annotation(title = sprintf("Visual checks — Study A  cond%03d rep%03d", cid, rid))

      ggsave(file.path(out_dir_cond, sprintf("checks_A_cond%03d_rep%03d.png", cid, rid)),
             g, width = 12, height = 14, dpi = 150)

    } else {
      z  <- sim$truth$innovations[-1, , drop = FALSE]
      u1 <- exp_cdf_dir(z[,1], dir_sign[1]); u2 <- exp_cdf_dir(z[,2], dir_sign[2])

      g <- ((p_hist(z[,1], "B: innovations ζ1", sprintf("cond%03d rep%03d", cid, rid)) |
             p_hist(z[,2], "B: innovations ζ2", sprintf("cond%03d rep%03d", cid, rid))) /
            (p_qq(z[,1], dir_sign[1], "QQ ζ1 vs Exp±") |
             p_qq(z[,2], dir_sign[2], "QQ ζ2 vs Exp±")) /
            (p_pit(u1, "PIT(u1)") | p_pit(u2, "PIT(u2)")) /
             p_cop(u1, u2, rho, "Gaussian-score scatter (innovation layer)")) +
        plot_annotation(title = sprintf("Visual checks — Study B  cond%03d rep%03d", cid, rid))

      ggsave(file.path(out_dir_cond, sprintf("checks_B_cond%03d_rep%03d.png", cid, rid)),
             g, width = 12, height = 14, dpi = 150)
    }
  }
  message("✓ Visual checks saved under: ", OUT_DIR)
}

run_post_sim_checks_sem <- function(DATA_DIR = "data", CHECKS_DIR = "checks_sem") {
  run_visual_checks_sem(DATA_DIR = DATA_DIR, OUT_DIR = CHECKS_DIR, reps_to_plot = 3)
}
