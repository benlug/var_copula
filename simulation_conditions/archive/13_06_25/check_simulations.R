###########################################################################
# check_simulations.R  – Visual diagnostics (ASCII‑only labels, no warnings)
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

reshape_phi <- function(arr) {
  tibble(
    i = rep(seq_len(dim(arr)[1]), each = 4),
    element = rep(c("phi11", "phi12", "phi21", "phi22"),
      times = dim(arr)[1]
    ),
    value = c(arr[, 1, 1], arr[, 1, 2], arr[, 2, 1], arr[, 2, 2])
  )
}

run_post_sim_checks_ml <- function(data_dir, checks_dir,
                                   n_subjects_to_plot = 3,
                                   n_reps_to_plot = 1) {
  files <- list.files(data_dir, "^sim_data_ml_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(files)) {
    message("No data to check.")
    return(invisible())
  }

  cond_ids <- sort(unique(sub(
    "^.*_cond(\\d+)_rep.*", "\\1",
    basename(files)
  )))

  for (cid in cond_ids) {
    pdf_path <- file.path(
      checks_dir,
      sprintf("checks_ml_cond_%03d.pdf", as.integer(cid))
    )
    pdf(pdf_path, 11, 8.5) # base PDF only knows ASCII fonts

    # ---------- title page (ASCII only) ---------------------------------
    grid.newpage()
    grid.text(sprintf("Diagnostics - condition %03d", as.integer(cid)),
      gp = gpar(fontsize = 18, fontface = "bold")
    )

    reps <- head(list.files(data_dir,
      sprintf("sim_data_ml_cond%03d_rep\\d+\\.rds", as.integer(cid)),
      full.names = TRUE
    ), n_reps_to_plot)

    for (rf in reps) {
      obj <- readRDS(rf)
      pick <- sample(
        seq_len(obj$N),
        min(n_subjects_to_plot, obj$N)
      )

      # ----- time series -------------------------------------------------
      p_ts <- obj$data |>
        filter(i %in% pick) |>
        pivot_longer(cols = c(y1, y2)) |>
        ggplot(aes(t, value, colour = name)) +
        geom_line() +
        facet_wrap(~i, ncol = 1) +
        theme_bw(9) +
        labs(
          title = sprintf(
            "Condition %03d  |  Rep %d",
            obj$condition_id, obj$rep_i
          ),
          x = "time", y = "value", colour = ""
        )
      print(p_ts)

      # ----- Phi histogram ----------------------------------------------
      p_phi <- reshape_phi(obj$phi_subject) |>
        ggplot(aes(value)) +
        geom_histogram(bins = 30) +
        facet_wrap(~element, scales = "free") +
        theme_bw(9) +
        labs(title = "Subject-specific Phi coefficients")
      print(p_phi)
    }

    dev.off()
    message("✓ Diagnostics written:", basename(pdf_path))
  }
}
