###############################################################################
# plot_divergences.R
#
# This script creates a histogram of the number of divergent transitions
# for each simulation condition, at the replication level.
# One page per condition in a single PDF file.
###############################################################################

# --- Libraries ---
library(tidyverse)
library(grid)
library(gridExtra)
library(forcats)
library(stringr)

# --- Configuration ---
# Adjust these paths as needed so they match your folder structure:
RESULTS_DIR <- getwd() # e.g., "path/to/results"
DATA_DIR <- file.path(RESULTS_DIR, "../data")
OUTPUT_FILE <- file.path(RESULTS_DIR, "divergence_plots_all_conditions.pdf")

# Filenames for condition data and sampler summaries
sim_conds_file <- file.path(DATA_DIR, "sim_conditions.rds")
sampler_info_file <- file.path(RESULTS_DIR, "sampler_summary.rds")

# --- 1. Load data ---
if (!file.exists(sim_conds_file)) {
  stop("Cannot find sim_conditions.rds at: ", sim_conds_file)
}
sim_conditions_df <- readRDS(sim_conds_file) %>%
  distinct(condition_id, .keep_all = TRUE) %>%
  mutate(condition_id = as.integer(condition_id))

if (!file.exists(sampler_info_file)) {
  stop("Cannot find sampler_summary.rds at: ", sampler_info_file)
}
sampler_info_df <- readRDS(sampler_info_file) %>%
  mutate(condition_id = as.integer(condition_id))

# --- 2. Merge basic condition info into sampler data (optional) ---
# This is just so we can label plots more clearly if needed.
sampler_info_df <- sampler_info_df %>%
  left_join(sim_conditions_df, by = "condition_id")

# --- 3. Create a PDF with one page per condition ---
pdf(OUTPUT_FILE, width = 8, height = 6)

all_condition_ids <- sort(unique(sampler_info_df$condition_id))

for (cond_id in all_condition_ids) {
  cond_data <- sampler_info_df %>%
    filter(condition_id == cond_id)

  # If no data for this condition, skip
  if (nrow(cond_data) == 0) next

  # Create a histogram of divergences across replications.
  # You could facet by model code if you want, or color them, etc.
  # The example below facets by fitted_model_code for clarity.
  p <- ggplot(cond_data, aes(x = divergences)) +
    geom_histogram(
      binwidth = 1,
      fill = "skyblue", color = "black", alpha = 0.7
    ) +
    facet_wrap(~fitted_model_code, ncol = 2, scales = "free_y") +
    labs(
      title = paste0("Number of Divergent Transitions\nCondition ", cond_id),
      subtitle = paste(
        "Histogram of divergences across all replications and models.",
        "One bin = 1 divergence. Zero is ideal."
      ),
      x = "Number of Divergences",
      y = "Count of Replications"
    ) +
    theme_bw(base_size = 10)

  # Print the plot (one page per condition)
  grid.newpage()
  print(p)
}

dev.off()

cat("Divergence plots saved to:", OUTPUT_FILE, "\n")
