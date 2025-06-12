###########################################################################
# analysis_singlelevel.R
#
# • Assumes 192×100 datasets have been simulated and fitted
# • Reads fits/, computes parameter bias, SD‑Bias, coverage
# • Outputs CSVs + quick convergence plots
###########################################################################

library(rstan)
library(tidyverse)
library(posterior)
library(bayesplot)

BASE_DIR <- this.dir() # folder with run_pipeline.R
DATA_DIR <- file.path(BASE_DIR, "data")
FITS_DIR <- file.path(BASE_DIR, "fits")
RESULTS_DIR <- file.path(BASE_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE)
QC_DIR <- file.path(RESULTS_DIR, "quick_checks")
dir.create(QC_DIR, showWarnings = FALSE)

design <- readRDS(file.path(DATA_DIR, "sim_conditions_singlelevel.rds")) %>%
  select(condition_id, skew_level, direction, T, rho, VARset)

# -------------------------------------------------------------------------
# helper: true parameter lookup
# -------------------------------------------------------------------------
true_value <- function(name, sim_obj) {
  if (name == "rho") {
    return(sim_obj$rho)
  }
  if (grepl("^phi", name)) {
    idx <- as.integer(substr(name, 4, 4))
    jdx <- as.integer(substr(name, 5, 5))
    return(sim_obj$phi_matrix[idx, jdx])
  }
  if (name == "mu[1]") {
    return(0)
  }
  if (name == "mu[2]") {
    return(0)
  }
  if (name == "sigma[1]" || name == "sigma[2]") {
    return(1)
  } # Var≈1 by design
  NA_real_
}

# valid parameters of interest
keep_params <- c(
  "mu[1]", "mu[2]",
  "phi11", "phi12", "phi21", "phi22",
  "sigma[1]", "sigma[2]",
  "rho"
)

# -------------------------------------------------------------------------
# iterate over fit files
# -------------------------------------------------------------------------
fit_files <- list.files(FITS_DIR,
  "^fit_(NG|SG)_cond\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)
stopifnot(length(fit_files) > 0)

param_res <- list()
sampler_res <- list()
pb <- utils::txtProgressBar(0, length(fit_files), style = 3)

for (i in seq_along(fit_files)) {
  fn <- fit_files[i]
  mat <- stringr::str_match(
    basename(fn),
    "^fit_([A-Z]{2})_cond(\\d+)_rep(\\d+)\\.rds$"
  )
  code <- mat[1, 2]
  cond <- as.integer(mat[1, 3])
  rep <- as.integer(mat[1, 4])

  fit <- readRDS(fn)
  sim_obj <- readRDS(
    file.path(DATA_DIR, sprintf("sim_data_cond%03d_rep%03d.rds", cond, rep))
  )

  s <- summary(fit)$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    filter(parameter %in% keep_params) %>%
    transmute(
      parameter,
      post_mean = mean,
      post_sd   = sd,
      ci_low    = `2.5%`,
      ci_high   = `97.5%`
    )
  s$true_val <- map_dbl(s$parameter, true_value, sim_obj = sim_obj)
  s$bias <- s$post_mean - s$true_val
  s$rel_bias <- ifelse(abs(s$true_val) > 1e-8,
    s$bias / abs(s$true_val), NA_real_
  )
  s$coverage <- s$true_val >= s$ci_low & s$true_val <= s$ci_high
  s$condition_id <- cond
  s$rep_i <- rep
  s$fitted_model <- code
  param_res[[length(param_res) + 1]] <- s

  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  div <- sum(map_dbl(sp, ~ sum(.x[, "divergent__"])))
  eF <- mean(map_dbl(sp, ~ {
    e <- .x[, "energy__"]
    n <- length(e)
    if (n < 3) {
      return(NA_real_)
    }
    em <- mean(e)
    varE <- var(e)
    ac1 <- sum((e[-1] - em) * (e[-n] - em)) / (n - 1)
    varE / (varE + 2 * ac1)
  }), na.rm = TRUE)
  sampler_res[[length(sampler_res) + 1]] <- tibble(
    condition_id = cond, rep_i = rep, model = code,
    divergences = div, eFMI = eF
  )
  utils::setTxtProgressBar(pb, i)
}
close(pb)

param_df <- bind_rows(param_res)
sampler_df <- bind_rows(sampler_res)

# ------------------------------ SD‑Bias per condition ---------------------
sd_bias_df <- param_df %>%
  group_by(condition_id, fitted_model, parameter) %>%
  summarise(
    emp_sd = sd(post_mean),
    mean_psd = mean(post_sd),
    sd_bias = mean_psd - emp_sd,
    .groups = "drop"
  )

# merge design factors
param_df <- param_df %>% left_join(design, by = "condition_id")
sampler_df <- sampler_df %>% left_join(design, by = "condition_id")
sd_bias_df <- sd_bias_df %>% left_join(design, by = "condition_id")

# ------------------------------ write CSVs --------------------------------
write_csv(
  param_df,
  file.path(RESULTS_DIR, "parameter_summary.csv")
)
write_csv(
  sampler_df,
  file.path(RESULTS_DIR, "sampler_summary.csv")
)
write_csv(
  sd_bias_df,
  file.path(RESULTS_DIR, "sd_bias_summary.csv")
)
cat("CSV files written in results/\n")

# ------------------------------ quick PDFs --------------------------------
theme_set(theme_bw(base_size = 9))
for (cid in sort(unique(sampler_df$condition_id))) {
  pdf(file.path(QC_DIR, sprintf("quick_cond_%03d.pdf", cid)),
    width = 7, height = 4
  )
  d <- sampler_df %>% filter(condition_id == cid)
  if (nrow(d) == 0) {
    dev.off()
    next
  }
  p1 <- ggplot(d, aes(model, divergences, fill = model)) +
    geom_violin(scale = "width", alpha = .6) +
    labs(
      title = sprintf("Condition %03d – divergences", cid),
      x = NULL, y = "# divergent transitions"
    ) +
    scale_y_sqrt() +
    guides(fill = "none")
  p2 <- ggplot(d, aes(model, eFMI, fill = model)) +
    geom_violin(scale = "width", alpha = .6) +
    geom_hline(yintercept = 0.3, linetype = "dashed", colour = "red") +
    labs(title = "E‑FMI", x = NULL, y = "Energy fraction MI") +
    guides(fill = "none")
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  dev.off()
}
cat("Quick convergence PDFs in results/quick_checks/\n")
