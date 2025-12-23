###########################################################################
# check_simulations.R
# Visual sanity checks - produces one PDF per condition (max 5 reps each)
# Improved layout with PIT diagnostics and copula correlation verification
###########################################################################

# Check for required packages
if (!requireNamespace("sn", quietly = TRUE)) {
  message("Package 'sn' not available - PIT plots will use empirical CDF")
}

## ========================================================================
## Helper Functions
## ========================================================================

# Compute residuals using TRUE generating parameters
compute_residuals_TRUE <- function(y, mu, phi) {
  Tn <- nrow(y)
  if (Tn < 2) return(list(res1 = numeric(0), res2 = numeric(0)))
  
  mu_matrix <- matrix(mu, nrow = Tn - 1, ncol = 2, byrow = TRUE)
  predicted <- mu_matrix + y[1:(Tn - 1), ] %*% t(phi)
  res <- y[2:Tn, ] - predicted
  list(res1 = res[, 1], res2 = res[, 2])
}

# Compute PIT values for residuals
compute_pit <- function(residuals, margin_info) {
  if (is.null(margin_info) || length(residuals) == 0) {
    return(NULL)
  }
  
  if (margin_info$type == "skewnormal") {
    if (!requireNamespace("sn", quietly = TRUE)) {
      # Fallback to empirical CDF
      return(rank(residuals) / (length(residuals) + 1))
    }
    alpha <- margin_info$alpha
    delta <- alpha / sqrt(1 + alpha^2)
    omega <- sqrt(1 / (1 - 2 * delta^2 / pi))
    xi <- -omega * delta * sqrt(2 / pi)
    return(sn::psn(residuals, xi = xi, omega = omega, alpha = alpha))
    
  } else if (margin_info$type == "chisq") {
    # For chi-squared: reverse the standardization
    df <- margin_info$df
    mean_chi <- df
    sd_chi <- sqrt(2 * df)
    mirror <- isTRUE(margin_info$mirror)
    
    if (mirror) {
      # If mirrored: residuals = -(x_std), and x_std = (x_raw - mean)/sd
      # So x_raw = mean - residuals * sd
      # But we also used 1-u in generation, so PIT = 1 - pchisq(x_raw, df)
      x_raw <- mean_chi - residuals * sd_chi
      # Clamp to valid range for pchisq
      x_raw <- pmax(0, x_raw)
      return(1 - pchisq(x_raw, df = df))
    } else {
      x_raw <- residuals * sd_chi + mean_chi
      x_raw <- pmax(0, x_raw)
      return(pchisq(x_raw, df = df))
    }
    
  } else {
    # Fallback to empirical CDF
    return(rank(residuals) / (length(residuals) + 1))
  }
}

# Color palette
pal <- list(
  series = c("#2166AC", "#B2182B"),
  hist = "#4393C3",
  pit = "#762A83",
  qq = "#1A1A1A",
  acf = "#4393C3",
  copula_ok = "#4DAF4A",
  copula_warn = "#E41A1C"
)

## ========================================================================
## Copula Correlation Verification
## ========================================================================

verify_residual_correlation <- function(res1, res2, expected_rho) {
  if (length(res1) < 10 || length(res2) < 10) {
    return(list(
      pearson = NA, kendall = NA, spearman = NA,
      tau_expected = NA, status = "insufficient_data"
    ))
  }
  
  # Empirical correlations
  r_pearson <- cor(res1, res2, method = "pearson")
  tau_kendall <- cor(res1, res2, method = "kendall")
  r_spearman <- cor(res1, res2, method = "spearman")
  
 # Expected Kendall's tau for Gaussian copula: tau = (2/pi) * arcsin(rho)
  tau_expected <- (2/pi) * asin(expected_rho)
  
  # Check sign consistency (most important)
  sign_ok <- sign(tau_kendall) == sign(expected_rho) || abs(expected_rho) < 0.1
  
  # Check magnitude (with tolerance for finite sample)
  n <- length(res1)
  tau_se <- sqrt(2 * (2*n + 5) / (9 * n * (n-1)))  # Approximate SE
  magnitude_ok <- abs(tau_kendall - tau_expected) < 3 * tau_se
  
  status <- if (sign_ok && magnitude_ok) "ok" else if (sign_ok) "magnitude_warning" else "sign_error"
  
  list(
    pearson = r_pearson,
    kendall = tau_kendall,
    spearman = r_spearman,
    tau_expected = tau_expected,
    tau_se = tau_se,
    status = status
  )
}

## ========================================================================
## Main Plotting Function
## ========================================================================

plot_diagnostics <- function(df, true_par, cond_lbl, rep_lbl) {
  ymat <- as.matrix(df[, c("y1", "y2")])
  res <- compute_residuals_TRUE(ymat, true_par$mu, true_par$phi)
  
  # Compute PITs
  pit1 <- compute_pit(res$res1, true_par$margin1)
  pit2 <- compute_pit(res$res2, true_par$margin2)
  
  # Verify copula correlation
  cop_check <- verify_residual_correlation(res$res1, res$res2, true_par$copula_rho)
  
  # Setup layout: 4 rows x 4 cols
  old <- par(mfrow = c(4, 4), 
             mar = c(3.5, 3.5, 2, 1), 
             mgp = c(2, 0.6, 0),
             oma = c(0, 0, 3, 0),
             cex.main = 0.9,
             cex.lab = 0.85,
             cex.axis = 0.8)
  on.exit(par(old))
  
  # ── Row 1: Time series and scatter ──
  # y1 time series
  plot(df$t, df$y1, type = "l", col = pal$series[1], lwd = 1.2,
       main = expression(paste("Time Series: ", Y[1])),
       xlab = "Time", ylab = expression(Y[1]))
  abline(h = 0, lty = 2, col = "gray60")
  
  # y2 time series
  plot(df$t, df$y2, type = "l", col = pal$series[2], lwd = 1.2,
       main = expression(paste("Time Series: ", Y[2])),
       xlab = "Time", ylab = expression(Y[2]))
  abline(h = 0, lty = 2, col = "gray60")
  
  # Scatter plot of observed
  plot(df$y1, df$y2, pch = 19, col = adjustcolor("black", 0.4), cex = 0.8,
       main = expression(paste("Scatter: ", Y[1], " vs ", Y[2])),
       xlab = expression(Y[1]), ylab = expression(Y[2]))
  if (nrow(df) > 2) {
    r <- cor(df$y1, df$y2)
    legend("topleft", bty = "n", cex = 0.8,
           legend = sprintf("r = %.2f", r))
  }
  
  # Scatter of residuals with copula check
  if (length(res$res1) > 0) {
    status_col <- switch(cop_check$status,
                         "ok" = pal$copula_ok,
                         "magnitude_warning" = "orange",
                         "sign_error" = pal$copula_warn,
                         "gray50")
    
    plot(res$res1, res$res2, pch = 19, col = adjustcolor("black", 0.4), cex = 0.8,
         main = expression(paste("Residuals: ", epsilon[1], " vs ", epsilon[2])),
         xlab = expression(epsilon[1]), ylab = expression(epsilon[2]))
    
    # Add correlation info with status color
    legend("topleft", bty = "n", cex = 0.75,
           text.col = c("black", status_col),
           legend = c(sprintf("tau = %.2f (exp: %.2f)", 
                              cop_check$kendall, cop_check$tau_expected),
                      toupper(cop_check$status)))
  } else {
    plot.new()
  }
  
  # ── Row 2: Residual histograms and QQ plots ──
  if (length(res$res1) > 0) {
    # Histogram res1
    hist(res$res1, breaks = "fd", col = adjustcolor(pal$hist, 0.6), border = "white",
         main = expression(paste("Histogram: ", epsilon[1])),
         xlab = expression(epsilon[1]), freq = FALSE)
    curve(dnorm(x), add = TRUE, col = "red", lwd = 1.5)
    
    # Histogram res2
    hist(res$res2, breaks = "fd", col = adjustcolor(pal$hist, 0.6), border = "white",
         main = expression(paste("Histogram: ", epsilon[2])),
         xlab = expression(epsilon[2]), freq = FALSE)
    curve(dnorm(x), add = TRUE, col = "red", lwd = 1.5)
    
    # QQ plot res1
    qqnorm(res$res1, main = expression(paste("Q-Q Normal: ", epsilon[1])),
           pch = 19, col = adjustcolor(pal$qq, 0.5), cex = 0.7)
    qqline(res$res1, col = "red", lwd = 1.5)
    
    # QQ plot res2
    qqnorm(res$res2, main = expression(paste("Q-Q Normal: ", epsilon[2])),
           pch = 19, col = adjustcolor(pal$qq, 0.5), cex = 0.7)
    qqline(res$res2, col = "red", lwd = 1.5)
  } else {
    for (j in 1:4) plot.new()
  }
  
  # ── Row 3: ACF plots and PIT histograms ──
  if (length(res$res1) > 0) {
    # ACF res1
    acf(res$res1, main = expression(paste("ACF: ", epsilon[1])),
        col = pal$acf, na.action = na.pass, ci.col = "gray50")
    
    # ACF res2
    acf(res$res2, main = expression(paste("ACF: ", epsilon[2])),
        col = pal$acf, na.action = na.pass, ci.col = "gray50")
    
    # PIT histogram res1
    if (!is.null(pit1)) {
      hist(pit1, breaks = seq(0, 1, by = 0.1), 
           col = adjustcolor(pal$pit, 0.6), border = "white",
           main = expression(paste("PIT: ", epsilon[1])),
           xlab = "PIT", freq = FALSE, xlim = c(0, 1))
      abline(h = 1, col = "red", lwd = 1.5, lty = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "PIT unavailable", cex = 0.9)
    }
    
    # PIT histogram res2
    if (!is.null(pit2)) {
      hist(pit2, breaks = seq(0, 1, by = 0.1),
           col = adjustcolor(pal$pit, 0.6), border = "white",
           main = expression(paste("PIT: ", epsilon[2])),
           xlab = "PIT", freq = FALSE, xlim = c(0, 1))
      abline(h = 1, col = "red", lwd = 1.5, lty = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "PIT unavailable", cex = 0.9)
    }
  } else {
    for (j in 1:4) plot.new()
  }
  
  # ── Row 4: Summary statistics ──
  if (length(res$res1) > 0) {
    # res1 summary
    plot.new()
    stats1 <- c(
      sprintf("Mean: %.3f", mean(res$res1)),
      sprintf("SD: %.3f", sd(res$res1)),
      sprintf("Skew: %.3f", moments_skewness(res$res1)),
      sprintf("Kurt: %.3f", moments_kurtosis(res$res1))
    )
    text(0.5, 0.7, expression(bold(paste("Summary: ", epsilon[1]))), cex = 1)
    text(0.5, 0.45, paste(stats1, collapse = "\n"), cex = 0.85)
    
    # res2 summary
    plot.new()
    stats2 <- c(
      sprintf("Mean: %.3f", mean(res$res2)),
      sprintf("SD: %.3f", sd(res$res2)),
      sprintf("Skew: %.3f", moments_skewness(res$res2)),
      sprintf("Kurt: %.3f", moments_kurtosis(res$res2))
    )
    text(0.5, 0.7, expression(bold(paste("Summary: ", epsilon[2]))), cex = 1)
    text(0.5, 0.45, paste(stats2, collapse = "\n"), cex = 0.85)
    
    # PIT uniformity test
    plot.new()
    if (!is.null(pit1) && !is.null(pit2)) {
      ks1 <- ks.test(pit1, "punif")
      ks2 <- ks.test(pit2, "punif")
      text(0.5, 0.8, expression(bold("PIT Uniformity (KS)")), cex = 1)
      text(0.5, 0.55, sprintf("PIT1: D=%.3f, p=%.3f", ks1$statistic, ks1$p.value), cex = 0.85)
      text(0.5, 0.35, sprintf("PIT2: D=%.3f, p=%.3f", ks2$statistic, ks2$p.value), cex = 0.85)
    } else {
      text(0.5, 0.5, "KS test unavailable", cex = 0.9)
    }
    
    # Copula verification panel
    plot.new()
    text(0.5, 0.9, expression(bold("Copula Verification")), cex = 1)
    
    cop_info <- c(
      sprintf("True rho: %.2f", true_par$copula_rho),
      sprintf("Pearson r: %.3f", cop_check$pearson),
      sprintf("Kendall tau: %.3f", cop_check$kendall),
      sprintf("Expected tau: %.3f", cop_check$tau_expected)
    )
    
    status_col <- switch(cop_check$status,
                         "ok" = pal$copula_ok,
                         "magnitude_warning" = "orange", 
                         "sign_error" = pal$copula_warn,
                         "black")
    
    text(0.5, 0.55, paste(cop_info, collapse = "\n"), cex = 0.8)
    text(0.5, 0.15, sprintf("Status: %s", toupper(cop_check$status)), 
         col = status_col, font = 2, cex = 0.9)
  } else {
    for (j in 1:4) plot.new()
  }
  
  # Main title
  mtext(sprintf("%s  |  Rep %d", cond_lbl, rep_lbl),
        outer = TRUE, line = 1, cex = 1.2, font = 2)
}

# Simple moment calculations (avoid dependency on moments package)
moments_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)
  m <- mean(x)
  s <- sd(x)
  sum((x - m)^3) / (n * s^3)
}

moments_kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA)
  m <- mean(x)
  s <- sd(x)
  sum((x - m)^4) / (n * s^4) - 3  # excess kurtosis
}

## ========================================================================
## Main Function
## ========================================================================

run_post_sim_checks_var1 <- function(data_dir, checks_dir) {
  dir.create(checks_dir, FALSE, TRUE)
  
  files <- list.files(data_dir, "^sim_data_cond\\d+_rep\\d+\\.rds$",
                      full.names = TRUE)
  if (!length(files)) {
    message("No simulation files found in ", data_dir)
    return(invisible())
  }

  cond_ids <- sort(unique(sub("^.*_cond(\\d+)_rep.*", "\\1", basename(files))))

  message("Generating diagnostic plots for ", length(cond_ids), " conditions...")
  
  # Also create a copula verification summary
  copula_checks <- list()

  for (cid in cond_ids) {
    cid_num <- as.integer(cid)
    fname_pdf <- file.path(checks_dir, sprintf("checks_cond_%03d.pdf", cid_num))

    # Find all replications for this condition
    reps <- list.files(data_dir,
                       sprintf("sim_data_cond%03d_rep\\d+\\.rds", cid_num),
                       full.names = TRUE)

    # Sample up to 5 replications
    reps_to_check <- sample(reps, min(5, length(reps)))

    if (length(reps_to_check) > 0) {
      pdf(fname_pdf, width = 14, height = 12)

      for (rf in reps_to_check) {
        dat <- try(readRDS(rf), silent = TRUE)
        if (inherits(dat, "try-error")) {
          message("Error reading: ", rf)
          next
        }

        # Create condition label with more info
        margin_type <- dat$true_params$margin1$type
        if (margin_type == "skewnormal") {
          alpha1 <- dat$true_params$margin1$alpha
          alpha2 <- dat$true_params$margin2$alpha
          margin_desc <- sprintf("SN(a=%+d,%+d)", alpha1, alpha2)
        } else {
          mirror1 <- if (isTRUE(dat$true_params$margin1$mirror)) "-" else "+"
          mirror2 <- if (isTRUE(dat$true_params$margin2$mirror)) "-" else "+"
          margin_desc <- sprintf("Chi2(%s,%s)", mirror1, mirror2)
        }
        
        cond_lbl <- sprintf("Cond %03d: %s, T=%d, rho=%.2f", 
                            cid_num, margin_desc, dat$T, dat$rho)

        plot_diagnostics(dat$data, dat$true_params, cond_lbl, dat$rep_i)
        
        # Store copula check for summary
        ymat <- as.matrix(dat$data[, c("y1", "y2")])
        res <- compute_residuals_TRUE(ymat, dat$true_params$mu, dat$true_params$phi)
        cop_check <- verify_residual_correlation(res$res1, res$res2, dat$true_params$copula_rho)
        
        copula_checks[[length(copula_checks) + 1]] <- data.frame(
          condition_id = cid_num,
          rep_id = dat$rep_i,
          margin_type = margin_type,
          true_rho = dat$true_params$copula_rho,
          tau_empirical = cop_check$kendall,
          tau_expected = cop_check$tau_expected,
          status = cop_check$status,
          stringsAsFactors = FALSE
        )
      }
      dev.off()
      message("Written: ", basename(fname_pdf))
    }
  }
  
  # Write copula verification summary
  if (length(copula_checks) > 0) {
    copula_summary <- do.call(rbind, copula_checks)
    write.csv(copula_summary, 
              file.path(checks_dir, "copula_verification_summary.csv"),
              row.names = FALSE)
    
    # Print summary
    message("\n=== Copula Verification Summary ===")
    status_table <- table(copula_summary$status)
    print(status_table)
    
    if (any(copula_summary$status == "sign_error")) {
      warning("SIGN ERRORS detected in copula correlation! Check simulation code.")
      sign_errors <- copula_summary[copula_summary$status == "sign_error", ]
      message("Conditions with sign errors:")
      print(unique(sign_errors[, c("condition_id", "true_rho", "tau_empirical")]))
    }
  }
}
