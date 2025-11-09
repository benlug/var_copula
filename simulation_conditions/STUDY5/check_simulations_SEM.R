###########################################################################
# check_simulations_SEM.R — Visual sanity checks for SEM skewness study
# Study B: checks state innovations (true layer).
# Study A: checks measurement residuals via EI fit if available; else y fallback.
# Saves PNGs to CHECKS_DIR/condXXX/.
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(patchwork)
  library(rstan)
})

# ---- helpers ------------------------------------------------------------
exp_cdf_dir <- function(e, s = 1, dir = +1L) {
  if (dir == +1L) out <- ifelse(e <= -s, 0, 1 - exp(-(e / s + 1)))
  else            out <- ifelse(e >=  s, 1, exp(e / s - 1))
  pmin(pmax(out, 1e-12), 1 - 1e-12)
}
nscore <- function(u) stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))

theme_sm <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10))

p_hist <- function(x, title, subtitle = "", bins = 40) {
  ggplot(data.frame(x=x), aes(x)) +
    geom_histogram(aes(y=after_stat(density)), bins=bins, alpha=.8) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "density") +
    theme_sm
}
p_pit <- function(u, title, subtitle = "", bins = 20) {
  ggplot(data.frame(u=u), aes(u)) +
    geom_histogram(bins=bins, color="white") +
    geom_hline(yintercept = length(u)/bins, linetype = 2) +
    labs(title = title, subtitle = subtitle, x = "PIT(u)", y = "count") +
    theme_sm
}
p_qq_exp_dir <- function(e, dir, s = 1, title = "QQ vs Exp±") {
  if (dir == -1L) e <- -e
  p  <- (seq_along(e) - 0.5)/length(e)
  qth <- -1 - log(1 - p)
  df  <- data.frame(theoretical = qth, sample = sort((e/s)))
  ggplot(df, aes(theoretical, sample)) +
    geom_abline(slope=1, intercept=0, linetype=2) +
    geom_point(alpha=.6, size=1.5) +
    labs(title = title, subtitle = ifelse(dir==1,"right-skew (+)","left-skew (−)"),
         x = "Theoretical (std Exp+)", y = "Sample (std & mirrored if needed)") +
    theme_sm
}
p_scatter_nscore <- function(u1, u2, rho, title = "Copula Gaussian scores") {
  df <- data.frame(z1 = nscore(u1), z2 = nscore(u2))
  cc <- suppressWarnings(cor(df$z1, df$z2, use = "complete.obs"))
  ggplot(df, aes(z1, z2)) +
    geom_point(alpha=.35, size=1) +
    labs(title = title,
         subtitle = sprintf("Empirical corr(z)=%.2f  (true rho=%.2f)", cc, rho),
         x = "Φ^{-1}(u1)", y = "Φ^{-1}(u2)") +
    theme_sm
}

# ---- main entry ---------------------------------------------------------
run_visual_checks_sem <- function(DATA_DIR = "data",
                                  FITS_DIR = "fits_sem",
                                  OUT_DIR  = "checks_sem",
                                  reps_to_plot = 3) {

  # Path safety: ensure character, non-empty, and create OUT_DIR
  if (missing(OUT_DIR) || is.null(OUT_DIR) || !is.character(OUT_DIR) || !nzchar(OUT_DIR[1])) {
    OUT_DIR <- "checks_sem"
  }
  if (missing(DATA_DIR) || is.null(DATA_DIR) || !is.character(DATA_DIR) || !nzchar(DATA_DIR[1])) {
    DATA_DIR <- "data"
  }
  if (missing(FITS_DIR) || is.null(FITS_DIR) || !is.character(FITS_DIR) || !nzchar(FITS_DIR[1])) {
    FITS_DIR <- "fits_sem"
  }
  OUT_DIR <- normalizePath(OUT_DIR, mustWork = FALSE)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  message("Visual-check output root → ", OUT_DIR)

  # locate simulations
  paths <- list.files(DATA_DIR, "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
  if (!length(paths)) {
    message("No simulation files found in ", DATA_DIR)
    return(invisible(NULL))
  }

  # parse ids robustly, discard non-matching files
  meta <- stringr::str_match(basename(paths), "cond(\\d+)_rep(\\d+)")
  df <- data.frame(
    path = paths,
    cid  = suppressWarnings(as.integer(meta[, 2])),
    rid  = suppressWarnings(as.integer(meta[, 3]))
  ) |>
    dplyr::filter(!is.na(cid), !is.na(rid)) |>
    dplyr::arrange(cid, rid) |>
    dplyr::group_by(cid) |>
    dplyr::slice_head(n = reps_to_plot) |>
    dplyr::ungroup()

  if (nrow(df) == 0) {
    message("No parsable sim files under ", DATA_DIR,
            " (expected: sim_data_cond###_rep###.rds)")
    return(invisible(NULL))
  }

  # draw per selected file
  for (k in seq_len(nrow(df))) {
    p   <- df$path[k]
    cid <- df$cid[k]; rid <- df$rid[k]

    out_dir_cond <- file.path(OUT_DIR, sprintf("cond%03d", as.integer(cid)))
    if (is.na(out_dir_cond) || !nzchar(out_dir_cond)) {
      message("Skipping file with unparsable condition id: ", p)
      next
    }
    dir.create(out_dir_cond, showWarnings = FALSE, recursive = TRUE)

    sim  <- readRDS(p)
    y    <- as.matrix(sim$data[,c("y1","y2")])
    Tt   <- nrow(y)
    rho  <- sim$rho
    dir_sign <- as.integer(sim$skew_signs)
    B    <- sim$phi_matrix
    study <- sim$sem_study

    if (identical(study, "B_latent")) {
      # state-layer checks
      ylag <- rbind(c(0,0), y[-Tt,,drop=FALSE])
      z1 <- y[,1] - (B[1,1]*ylag[,1] + B[1,2]*ylag[,2])
      z2 <- y[,2] - (B[2,1]*ylag[,1] + B[2,2]*ylag[,2])
      u1 <- exp_cdf_dir(z1, s = 1, dir = dir_sign[1])
      u2 <- exp_cdf_dir(z2, s = 1, dir = dir_sign[2])

      g1 <- p_hist(z1, "Study B: state innovations z1", sprintf("cond%03d rep%03d", cid, rid))
      g2 <- p_hist(z2, "Study B: state innovations z2", sprintf("cond%03d rep%03d", cid, rid))
      g3 <- p_qq_exp_dir(z1, dir_sign[1], 1, "QQ z1 vs Exp±")
      g4 <- p_qq_exp_dir(z2, dir_sign[2], 1, "QQ z2 vs Exp±")
      g5 <- p_pit(u1, "PIT(u1) for state innovations")
      g6 <- p_pit(u2, "PIT(u2) for state innovations")
      g7 <- p_scatter_nscore(u1, u2, rho, "Gaussian scores (state layer)")

      ((g1|g2)/(g3|g4)/(g5|g6)/g7) +
        plot_annotation(title = sprintf("Visual checks – Study B  cond%03d rep%03d", cid, rid)) -> gp
      ggsave(file.path(out_dir_cond, sprintf("checks_B_cond%03d_rep%03d.png", cid, rid)),
             gp, width = 12, height = 14, dpi = 150)

    } else if (identical(study, "A_indicator")) {
      # try EI fit
      fit_path <- file.path(FITS_DIR, sprintf("fit_EI_cond%03d_rep%03d.rds", cid, rid))
      fit <- NULL
      if (file.exists(fit_path)) {
        tmp <- readRDS(fit_path)
        if (inherits(tmp, "stanfit")) fit <- tmp
      }

      if (!is.null(fit)) {
        sm <- summary(fit, pars = c("state","eta","rho"))$summary
        state_mean <- matrix(NA_real_, nrow = Tt, ncol = 2)
        for (t in 1:Tt) for (j in 1:2) {
          rn <- sprintf("state[%d,%d]", t, j)
          state_mean[t,j] <- sm[rn, "mean"]
        }
        eta_mean  <- c(sm["eta[1]","mean"], sm["eta[2]","mean"])
        sigma_hat <- exp(eta_mean)
        rho_hat   <- if ("rho" %in% rownames(sm)) sm["rho","mean"] else NA_real_

        e1 <- y[,1] - state_mean[,1]
        e2 <- y[,2] - state_mean[,2]
        u1 <- exp_cdf_dir(e1, s = sigma_hat[1], dir = dir_sign[1])
        u2 <- exp_cdf_dir(e2, s = sigma_hat[2], dir = dir_sign[2])

        g1 <- p_hist(e1, "Study A: measurement residuals ê1", sprintf("cond%03d rep%03d", cid, rid))
        g2 <- p_hist(e2, "Study A: measurement residuals ê2", sprintf("cond%03d rep%03d", cid, rid))
        g3 <- p_qq_exp_dir(e1/sigma_hat[1], dir_sign[1], 1, "QQ ê1/s vs Exp±")
        g4 <- p_qq_exp_dir(e2/sigma_hat[2], dir_sign[2], 1, "QQ ê2/s vs Exp±")
        g5 <- p_pit(u1, "PIT(u1) for measurement residuals")
        g6 <- p_pit(u2, "PIT(u2) for measurement residuals")
        g7 <- p_scatter_nscore(u1, u2, rho,
                               sprintf("Gaussian scores (meas. layer)  (true=%.2f; fit≈%.2f)", rho, rho_hat))

        ((g1|g2)/(g3|g4)/(g5|g6)/g7) +
          plot_annotation(title = sprintf("Visual checks – Study A  cond%03d rep%03d", cid, rid)) -> gp
        ggsave(file.path(out_dir_cond, sprintf("checks_A_cond%03d_rep%03d.png", cid, rid)),
               gp, width = 12, height = 14, dpi = 150)

      } else {
        # fallback: y marginals + scatter
        sdir <- paste0(ifelse(dir_sign[1]==1,"+","-"), ifelse(dir_sign[2]==1,"+","-"))
        g1 <- p_hist(y[,1], "Study A fallback: y1", sprintf("dir=%s  cond%03d rep%03d", sdir, cid, rid))
        g2 <- p_hist(y[,2], "Study A fallback: y2", sprintf("dir=%s  cond%03d rep%03d", sdir, cid, rid))
        g3 <- ggplot(data.frame(y1=y[,1], y2=y[,2]), aes(y1,y2)) +
          geom_point(alpha=.35, size=1) +
          labs(title="Scatter(y1,y2)", subtitle="measurement-layer skew expected", x="y1", y="y2") + theme_sm
        ((g1|g2)/g3) +
          plot_annotation(title = sprintf("Visual checks – Study A (fallback)  cond%03d rep%03d", cid, rid)) -> gp
        ggsave(file.path(out_dir_cond, sprintf("checks_A_fallback_cond%03d_rep%03d.png", cid, rid)),
               gp, width = 12, height = 10, dpi = 150)
      }
    }
  }

  message("✓ Visual checks saved under: ", OUT_DIR)
  invisible(TRUE)
}

# Pipeline-facing wrapper
run_post_sim_checks_sem <- function(DATA_DIR = "data", FITS_DIR = "fits_sem", CHECKS_DIR = "checks_sem") {
  run_visual_checks_sem(DATA_DIR = DATA_DIR, FITS_DIR = FITS_DIR, OUT_DIR = CHECKS_DIR, reps_to_plot = 3)
}
