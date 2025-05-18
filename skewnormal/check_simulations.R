###########################################################################
# check_simulations.R (Generalized Version - Corrected Filtering)
#
# Performs visual diagnostic checks on simulated multilevel VAR(1) data.
# Generates multi-panel plots per subject for each replication/condition.
# Creates one PDF per condition containing plots for selected subjects.
# Called by run_pipeline.R
###########################################################################

# --- Libraries ---
# Although run_pipeline.R might load these, it's good practice to ensure
# they are available or checked for within the script/functions that use them.
# library(forecast) # Checked for within plot_checks_for_subject
# library(graphics) # Checked for within plot_checks_for_subject
# library(stats)    # Base package, usually available
# library(dplyr)    # Checked for within run_post_sim_checks
# library(grDevices)# For pdf() and dev.off() - base package


# --- Helper Functions ---

#' Compute Approximate Residuals using Base AR Parameters Only
#'
#' Calculates residuals assuming a simple VAR(1) model using only the
#' global base AR parameters, ignoring subject-specific means or deviations.
#' This provides a basic check on the simulated dynamics.
#'
#' @param y Matrix of bivariate time series data (T x 2).
#' @param phi_base Numeric vector of length 4: c(phi11, phi12, phi21, phi22).
#' @return A list containing two numeric vectors: res1 and res2 (residuals for y1 and y2).
#'
compute_residuals_baseAR <- function(y, phi_base) {
    if (!is.numeric(phi_base) || length(phi_base) != 4) {
        stop("phi_base must be a numeric vector of length 4.")
    }
    phi11_base <- phi_base[1]; phi12_base <- phi_base[2]
    phi21_base <- phi_base[3]; phi22_base <- phi_base[4]

    if (!is.matrix(y)) y <- as.matrix(y)
    if (ncol(y) != 2) stop("Input matrix 'y' must have 2 columns.")

    y1_vals <- y[, 1]; y2_vals <- y[, 2]
    T_val <- nrow(y)

    if (T_val < 2) {
        warning("Time series too short (T < 2) to calculate residuals.")
        return(list(res1 = numeric(0), res2 = numeric(0)))
    }

    res1 <- numeric(T_val - 1)
    res2 <- numeric(T_val - 1)

    for (t in 2:T_val) {
        # Using simple mean=0 assumption for base AR residuals here.
        mu1_pred <- phi11_base * y1_vals[t - 1] + phi12_base * y2_vals[t - 1]
        mu2_pred <- phi21_base * y1_vals[t - 1] + phi22_base * y2_vals[t - 1]
        res1[t - 1] <- y1_vals[t] - mu1_pred
        res2[t - 1] <- y2_vals[t] - mu2_pred
    }
    list(res1 = res1, res2 = res2)
}


#' Plot Diagnostic Checks for a Single Subject
#'
#' Creates a multi-panel plot showing time series, scatter plots,
#' approximate residual diagnostics (histograms, QQ plots, time plots),
#' and ACF/PACF for the original series.
#'
#' @param df_sub Data frame containing data for a single subject (columns i, t, y1, y2).
#' @param subject_id The ID of the subject being plotted.
#' @param true_params List containing true simulation parameters, must include phi_base.
#' @param cond_id The condition ID for labeling.
#' @param rep_id The replication ID for labeling.
#' @return NULL. Plots to the active graphics device.
#'
plot_checks_for_subject <- function(df_sub, subject_id, true_params, cond_id, rep_id) {
    # Ensure necessary graphics packages are loaded
    if (!requireNamespace("graphics", quietly = TRUE)) stop("Package 'graphics' needed for plotting.")
    if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' needed for qqnorm, etc.")

    # Check for forecast package for ACF/PACF, fallback to base if missing
    acf_fun <- stats::acf; pacf_fun <- stats::pacf
    if (!requireNamespace("forecast", quietly = TRUE)) {
        warning("Package 'forecast' not found, using base acf/pacf for checks.", call. = FALSE)
    } else {
        acf_fun <- forecast::Acf; pacf_fun <- forecast::Pacf
    }

    # Helper function to plot an error message on a panel
    plot_error_message <- function(msg) {
        plot.new() # Start new plot panel
        # Use tryCatch for safety, title() can fail if device invalid
        tryCatch(title(main = msg, col.main="red", line = -1, cex.main=0.9), error=function(e) {})
        box() # Add frame to make it visible
    }

    # --- Data Validation ---
    if (!is.data.frame(df_sub) || nrow(df_sub) < 2) {
        warning(sprintf("Subject %d (Cond %s, Rep %s) insufficient data (rows=%d).",
                        subject_id, cond_id, rep_id, nrow(df_sub)), call. = FALSE)
        plot_error_message(sprintf("Subject %d \nInsufficient Data", subject_id))
        return()
    }
    if (!all(c("t", "y1", "y2") %in% names(df_sub))) {
         warning(sprintf("Subject %d (Cond %s, Rep %s) missing required columns (t, y1, y2).",
                        subject_id, cond_id, rep_id), call. = FALSE)
        plot_error_message(sprintf("Subject %d \nMissing Columns", subject_id))
        return()
    }
    y_mat <- as.matrix(df_sub[, c("y1", "y2")])

    # --- Residual Calculation ---
    if (is.null(true_params$phi_base) || length(true_params$phi_base) != 4) {
        warning(sprintf("phi_base missing or incorrect length in true_params for Cond %s, Rep %s.",
                        cond_id, rep_id), call. = FALSE)
        plot_error_message(sprintf("Subject %d \nMissing phi_base", subject_id))
        return()
    }
    phi_base <- true_params$phi_base
    res_list <- tryCatch(compute_residuals_baseAR(y_mat, phi_base),
                         error = function(e) {
                             warning("Error computing residuals: ", e$message, call. = FALSE); NULL
                         })

    if (is.null(res_list) || length(res_list$res1) == 0 || length(res_list$res2) == 0) {
        warning(sprintf("Residual calculation failed for Subject %d (Cond %s, Rep %s).",
                        subject_id, cond_id, rep_id), call. = FALSE)
        plot_error_message(sprintf("Subject %d \nResidual Calc Failed", subject_id))
        return()
    }
    res1 <- res_list$res1; res2 <- res_list$res2

    # --- Setup Plot Layout ---
    # Save current parameters and restore on exit
    old_par <- par(mfrow = c(4, 4), mar = c(3, 3.5, 2.5, 1), oma = c(0, 0, 3, 0), mgp = c(2, 0.7, 0))
    on.exit(par(old_par), add = TRUE) # Ensure graphical parameters are reset

    # --- Generate Plots ---
    # Row 1: Time Series
    tryCatch(plot(df_sub$t, y_mat[, 1], type = "l", col = "blue", main = "Time Series y1", xlab = "Time", ylab = "y1"), error=function(e) plot_error_message("Plot y1 Error"))
    tryCatch(plot(df_sub$t, y_mat[, 2], type = "l", col = "red", main = "Time Series y2", xlab = "Time", ylab = "y2"), error=function(e) plot_error_message("Plot y2 Error"))
    plot.new(); box(); plot.new(); box() # Spacers with visible box

    # Row 2: Scatter & Residual Histograms
    tryCatch(plot(y_mat[, 1], y_mat[, 2], main = "Scatter (y1 vs y2)", xlab = "y1", ylab = "y2", pch = 19, col = adjustcolor("purple", alpha.f = 0.6)), error=function(e) plot_error_message("Scatter Error"))
    tryCatch(hist(res1, breaks = "Sturges", main = "Hist Approx Resid (y1)", xlab = "BaseAR Resid", col = "lightblue"), error=function(e) plot_error_message("Hist R1 Error"))
    tryCatch(hist(res2, breaks = "Sturges", main = "Hist Approx Resid (y2)", xlab = "BaseAR Resid", col = "lightcoral"), error=function(e) plot_error_message("Hist R2 Error"))
    plot.new(); box()

    # Row 3: Residual QQ Plots & Time Plots
    tryCatch({qqnorm(res1, main = "QQ Approx Resid (y1)", pch = 19, cex = 0.7); qqline(res1, col = "red", lwd = 2)}, error=function(e) plot_error_message("QQ R1 Error"))
    tryCatch({qqnorm(res2, main = "QQ Approx Resid (y2)", pch = 19, cex = 0.7); qqline(res2, col = "blue", lwd = 2)}, error=function(e) plot_error_message("QQ R2 Error"))
    tryCatch(plot(seq_along(res1), res1, type = "l", main = "Time Plot Approx Resid (y1)", xlab = "Index (t-1)", ylab = "BaseAR Resid", col = "grey30"), error=function(e) plot_error_message("Plot R1 Error")); abline(h = 0, col = "red", lty = 2)
    tryCatch(plot(seq_along(res2), res2, type = "l", main = "Time Plot Approx Resid (y2)", xlab = "Index (t-1)", ylab = "BaseAR Resid", col = "grey30"), error=function(e) plot_error_message("Plot R2 Error")); abline(h = 0, col = "blue", lty = 2)

    # Row 4: ACF/PACF of Original Series
    # Use the selected acf/pacf functions (forecast or base)
    tryCatch(acf_fun(y_mat[, 1], lag.max = 20, main = "ACF of y1", na.action=na.pass), error=function(e) plot_error_message("ACF y1 Error"))
    tryCatch(pacf_fun(y_mat[, 1], lag.max = 20, main = "PACF of y1", na.action=na.pass), error=function(e) plot_error_message("PACF y1 Error"))
    tryCatch(acf_fun(y_mat[, 2], lag.max = 20, main = "ACF of y2", na.action=na.pass), error=function(e) plot_error_message("ACF y2 Error"))
    tryCatch(pacf_fun(y_mat[, 2], lag.max = 20, main = "PACF of y2", na.action=na.pass), error=function(e) plot_error_message("PACF y2 Error"))

    # --- Add Overall Title ---
    title_text <- sprintf("Subject %d | Condition %s | Replication %s", subject_id, cond_id, rep_id)
    # Use tryCatch for mtext as well, in case of issues with outer margins
    tryCatch(mtext(title_text, outer = TRUE, line = 1, cex = 1.2, font = 2), error=function(e) {})

    invisible(NULL) # Return nothing explicitly
}


# --- Main Function to Process Files ---

#' Run post-simulation visual checks for all conditions and replications.
#' Creates one PDF per condition in the `checks_dir`.
#'
#' @param data_dir Directory containing the simulation RDS files (e.g., "data_sn").
#' @param checks_dir Directory where the output PDF check files will be saved (e.g., "checks_sn").
#' @return Invisibly returns NULL. Saves PDF files to disk.
#'
run_post_sim_checks <- function(data_dir, checks_dir) {
    # Ensure dplyr is available for filter()
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Package 'dplyr' needed for run_post_sim_checks function.")
    }
    # Ensure grDevices is available for PDF output
    if (!requireNamespace("grDevices", quietly = TRUE)) {
         stop("Package 'grDevices' needed for PDF output.")
    }

    # --- Validate Directories ---
    if (!dir.exists(data_dir)) {
        stop("Data directory not found: ", data_dir)
    }
    if (!dir.exists(checks_dir)) {
        warning("Checks directory not found, creating: ", checks_dir, call. = FALSE)
        tryCatch(dir.create(checks_dir, recursive = TRUE),
                 error = function(e) stop("Failed to create checks directory: ", e$message))
    }

    # --- Find and Parse Data Files ---
    data_files <- list.files(data_dir, pattern = "^sim_data_cond\\d+_rep\\d+\\.rds$", full.names = TRUE)
    if (length(data_files) == 0) {
        cat(sprintf("No simulation data files found matching pattern in '%s'.\n", data_dir))
        return(invisible(NULL))
    }

    # Extract unique condition IDs from filenames (robustly handles padding)
    cond_ids_chr <- unique(sub("^.*_cond(\\d+)_rep\\d+\\.rds$", "\\1", basename(data_files)))
    cond_ids <- sort(as.integer(cond_ids_chr[!is.na(as.integer(cond_ids_chr))])) # Get unique numeric IDs

    if (length(cond_ids) == 0) {
        cat(sprintf("Could not extract valid numeric condition IDs from filenames in '%s'.\n", data_dir))
        return(invisible(NULL))
    }

    cat(sprintf("Found data for %d conditions in '%s': %s\n", length(cond_ids), data_dir, paste(cond_ids, collapse=", ")))

    # --- Process Each Condition ---
    for (cond_id in cond_ids) {
        cat(sprintf("\nProcessing Condition ID %d...\n", cond_id))

        # Filter files for the current condition using correct padding in the pattern
        cond_pattern <- sprintf("_cond%03d_rep", cond_id) # Assumes %03d padding
        cond_files <- data_files[grepl(cond_pattern, basename(data_files), fixed = FALSE)] # Use regex matching

        if (length(cond_files) == 0) {
            cat(sprintf("  Warning: No files found matching pattern '%s' for Condition ID %d.\n", cond_pattern, cond_id));
            next # Skip to the next condition
        }

        # --- Setup PDF Output for this Condition ---
        pdf_file <- file.path(checks_dir, sprintf("visual_checks_condition_%03d.pdf", cond_id))
        tryCatch({
            grDevices::pdf(pdf_file, width = 16, height = 10)
            pdf_opened <- TRUE
        }, error = function(e) {
            cat(sprintf("  ERROR: Failed to open PDF device for %s. Skipping condition. Error: %s\n", pdf_file, e$message))
            pdf_opened <- FALSE
        })
        if (!pdf_opened) next # Skip condition if PDF failed

        # --- Process Replications within the Condition ---
        cat(sprintf("  Processing %d file(s) for condition %d...\n", length(cond_files), cond_id))
        files_processed_count <- 0
        for (f in cond_files) {
            # Extract replication number using the padded condition ID
            rep_pattern <- sprintf("(?<=_cond%03d_rep)(\\d+)(?=\\.rds$)", cond_id)
            rep_id_match <- regmatches(basename(f), regexpr(rep_pattern, basename(f), perl=TRUE))
            if (length(rep_id_match) == 0) {
                 cat("    WARNING: Could not extract replication number from:", basename(f), "\n"); next
            }
            rep_id <- as.integer(rep_id_match)

            # Optional: Print progress less frequently if many reps
            # if (files_processed_count %% 10 == 0) {
                 cat(sprintf("    Rep %d (File: %s)...\n", rep_id, basename(f)))
            # }

            # --- Load and Validate Data ---
            sim_dat <- NULL
            tryCatch({ sim_dat <- readRDS(f) }, error = function(e) {
                cat("      ERROR reading file:", e$message, "\n"); sim_dat <<- NULL
            })
            if (is.null(sim_dat)) { next }

            # Check essential components
            # --- Validate Data Structure ---
            valid_structure <- FALSE # Default to FALSE
            if (is.list(sim_dat) &&
                !is.null(sim_dat$data) && is.data.frame(sim_dat$data) &&
                !is.null(sim_dat$true_params) && is.list(sim_dat$true_params)) {
                # Check columns in the data frame
                if (all(c("i", "t", "y1", "y2") %in% names(sim_dat$data))) {
                    # Check specific true parameter element
                    if (!is.null(sim_dat$true_params$phi_base) && length(sim_dat$true_params$phi_base) == 4) {
                        valid_structure <- TRUE
                    } else {
                         cat("      INFO: Structure check failed (missing or wrong length true_params$phi_base).\n")
                    }
                } else {
                     cat("      INFO: Structure check failed (missing columns i, t, y1, or y2 in sim_dat$data).\n")
                }
            } else {
                 cat("      INFO: Structure check failed (sim_dat not list or missing data/true_params).\n")
            }

            if (!valid_structure) {
                cat("      WARNING: Skipping file due to missing/invalid content or structure.\n")
                next # Skip to the next file
            }
            # --- End Validation ---

            # --- Extract Data and Parameters ---
            true_params <- sim_dat$true_params
            all_data <- sim_dat$data
            sub_ids <- unique(all_data$i)
            N_subj <- length(sub_ids)

            # --- Plot Checks for a Subset of Subjects ---
            # Plot max 5 subjects to keep PDFs manageable
            n_subjects_to_plot <- min(N_subj, 5)
            subjects_to_plot <- if (N_subj > n_subjects_to_plot) sample(sub_ids, n_subjects_to_plot) else sub_ids

            # cat(sprintf("      Plotting checks for %d subjects: %s\n", n_subjects_to_plot, paste(subjects_to_plot, collapse=", ")))
            plot_success <- TRUE
            for (s_id in subjects_to_plot) {
                tryCatch({
                     df_sub <- dplyr::filter(all_data, i == s_id)
                     plot_checks_for_subject(
                        df_sub = df_sub, subject_id = s_id, true_params = true_params,
                        cond_id = cond_id, rep_id = rep_id
                     )
                }, error = function(e) {
                     cat(sprintf("      ERROR plotting subject %d for Rep %d: %s\n", s_id, rep_id, e$message))
                     plot_success <<- FALSE # Assign to outer scope
                })
                if(!plot_success) break # Stop plotting for this rep if one subject fails
            } # End subject loop

            files_processed_count <- files_processed_count + 1

        } # End replication loop for this condition

        # --- Close PDF Device ---
        tryCatch(grDevices::dev.off(), error = function(e){
            cat("      Warning: Error closing PDF device. It might be closed already.\n")
        })
        cat(sprintf("  Finished plots for condition %d => %s\n", cond_id, pdf_file))
    } # End condition loop

    cat(sprintf("\nDone with all simulation checks. Output PDFs are in '%s/' directory.\n", checks_dir))
    invisible(NULL) # Return nothing
}

# --- Example Call (if run standalone) ---
# If you were to run this script directly (not recommended for the pipeline):
# DATA_DIR <- "data_sn" # Set path to your data directory
# CHECKS_DIR <- "checks_sn" # Set path to your desired output directory
# run_post_sim_checks(data_dir = DATA_DIR, checks_dir = CHECKS_DIR)
