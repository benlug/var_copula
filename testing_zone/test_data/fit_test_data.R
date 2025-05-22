# fit_test_data.R
# Fits VAR(1) copula models to the example dataset in testing_zone/test_data
# Saves simple diagnostic plots and RDS objects with fit results

# --- Libraries ---
library(rstan)
library(dplyr)
library(ggplot2)

# Determine script directory; fall back to current working directory if
# this.path is unavailable (e.g., when running interactively)
SCRIPT_DIR <- tryCatch(this.path::this.dir(), error = function(e) getwd())
setwd(SCRIPT_DIR)
ROOT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
# --- load data ---
df <- read.table("Example_skew_time_series.txt", header = TRUE)

# remove rows with missing values
clean_df <- na.omit(df)

# basic distribution diagnostics
zero_prop_var1 <- mean(clean_df$VAR1 == 0)
zero_prop_var2 <- mean(clean_df$VAR2 == 0)
cat(sprintf("Proportion of zeros - VAR1: %.3f, VAR2: %.3f\n", zero_prop_var1, zero_prop_var2))

p_var1 <- ggplot(clean_df, aes(x = VAR1)) +
  geom_histogram(bins = 50, color = "black", fill = "skyblue") +
  ggtitle("Distribution of VAR1")

p_var2 <- ggplot(clean_df, aes(x = VAR2)) +
  geom_histogram(bins = 50, color = "black", fill = "salmon") +
  ggtitle("Distribution of VAR2")

# save plots to the directory one level up from this script
plot_dir <- ROOT_DIR
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

ggsave(file.path(plot_dir, "VAR1_distribution.png"), p_var1, width = 5, height = 4)

ggsave(file.path(plot_dir, "VAR2_distribution.png"), p_var2, width = 5, height = 4)

# example time series plot for first ID
sample_id <- unique(clean_df$ID)[1]
example_df <- dplyr::filter(clean_df, ID == sample_id)

p_ts1 <- ggplot(example_df, aes(x = seq_along(VAR1), y = VAR1)) +
  geom_line() +
  labs(x = "Time", y = "VAR1", title = paste("VAR1 time series - ID", sample_id))

ggsave(file.path(plot_dir, "VAR1_timeseries_id1.png"), p_ts1, width = 5, height = 4)

# --- compile stan models ---
stan_files <- list(
  NG = file.path(ROOT_DIR, "model_NG_sl.stan"),
  NC = file.path(ROOT_DIR, "model_NC_sl.stan"),
  SG = file.path(ROOT_DIR, "model_SNG_sl.stan"),
  SC = file.path(ROOT_DIR, "model_SNC_sl.stan")
)

models <- lapply(stan_files, stan_model)

# directory to save fit objects
fits_dir <- file.path(ROOT_DIR, "test_fits")
if (!dir.exists(fits_dir)) dir.create(fits_dir, recursive = TRUE)

# fit models separately for each ID
ids <- unique(clean_df$ID)
for (id in ids) {
  dat_id <- dplyr::filter(clean_df, ID == id)
  if (nrow(dat_id) < 2) next

  T_val <- nrow(dat_id)
  y_list <- lapply(1:T_val, function(t) c(dat_id$VAR1[t], dat_id$VAR2[t]))
  stan_data <- list(T = T_val, y = y_list)

  for (code in names(models)) {
    fit <- sampling(models[[code]],
      data = stan_data,
      iter = 2000, warmup = 1000, chains = 4, refresh = 0
    )
    out_file <- file.path(fits_dir, sprintf("fit_%s_id%03d.rds", code, id))
    saveRDS(fit, file = out_file)
  }
}
