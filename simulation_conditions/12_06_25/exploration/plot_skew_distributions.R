#!/usr/bin/env Rscript
# plot_skew_distributions.R  –  density + histogram of the three residual forms

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(sn) # skew‑normal generator
})

set.seed(2025)

# --- helper to obtain xi and omega so that Var = 1 -----------------------
sn_params <- function(alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(1 / (1 - 2 * delta^2 / pi)) # enforces Var = 1
  xi <- -omega * delta * sqrt(2 / pi)
  list(xi = xi, omega = omega)
}

# --- simulate 20 000 draws per distribution -----------------------------
draws <- list(
  "SN  alpha = +4" = rsn(20000,
    xi = sn_params(+4)$xi,
    omega = sn_params(+4)$omega,
    alpha = +4
  ),
  "SN  alpha = -4" = rsn(20000,
    xi = sn_params(-4)$xi,
    omega = sn_params(-4)$omega,
    alpha = -4
  ),
  "SN  alpha = +9" = rsn(20000,
    xi = sn_params(+9)$xi,
    omega = sn_params(+9)$omega,
    alpha = +9
  ),
  "SN  alpha = -9" = rsn(20000,
    xi = sn_params(-9)$xi,
    omega = sn_params(-9)$omega,
    alpha = -9
  ),
  "chi2 df=1 (right)" = rchisq(20000, df = 1),
  "chi2 df=1 (mirrored left)" = -rchisq(20000, df = 1)
)

df <- bind_rows(lapply(names(draws), function(nm) {
  tibble(value = draws[[nm]], dist = nm)
}))

# --- colour palette without symbols -------------------------------------
pal <- c(
  "SN  alpha = +4" = "#1b9e77",
  "SN  alpha = -4" = "#1b9e77",
  "SN  alpha = +9" = "#d95f02",
  "SN  alpha = -9" = "#d95f02",
  "chi2 df=1 (right)" = "#7570b3",
  "chi2 df=1 (mirrored left)" = "#7570b3"
)

p <- ggplot(df, aes(value)) +
  geom_histogram(
    aes(
      y = after_stat(density),
      fill = dist
    ),
    bins = 60, alpha = 0.25, colour = NA
  ) +
  geom_density(aes(colour = dist), linewidth = 0.8) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, name = "") +
  facet_wrap(~dist, scales = "free", ncol = 2) +
  theme_bw(base_size = 9) +
  labs(
    title = "Marginal residual distributions used in the simulation",
    x = "value", y = "density"
  )

ggsave("skew_distributions.pdf", p, width = 8, height = 6)
