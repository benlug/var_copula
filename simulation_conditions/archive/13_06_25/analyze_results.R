###########################################################################
# analyze_results.R  –  Summarise bias / coverage etc.
###########################################################################

analyze_results_ml <- function(data_dir, fits_dir, res_dir) {
  library(rstan)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  count_div <- function(f) {
    sum(vapply(
      rstan::get_sampler_params(f, FALSE),
      \(x) sum(x[, "divergent__"]), 0
    ))
  }

  grab <- function(f, p) {
    s <- summary(f, pars = p)$summary
    tibble(
      param = rownames(s), mean = s[, "mean"], sd = s[, "sd"],
      l95 = s[, "2.5%"], u95 = s[, "97.5%"]
    )
  }

  CORE <- c("mu0[1]", "mu0[2]", paste0("Phi0[", 1:4, "]"), "rho")
  EXTRA <- c("omega[1]", "omega[2]", "alpha[1]", "alpha[2]")

  fits <- list.files(fits_dir, "^fit_(SNG|NG)_cond\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  stopifnot(length(fits) > 0)

  rep_tbl <- purrr::map_dfr(fits, \(f){
    m <- str_match(basename(f), "fit_(SNG|NG)_cond(\\d+)_rep(\\d+)\\.rds")
    mod <- m[2]
    cid <- as.integer(m[3])
    rid <- as.integer(m[4])

    sim <- readRDS(file.path(
      data_dir,
      sprintf("sim_data_ml_cond%03d_rep%03d.rds", cid, rid)
    ))
    truth <- list(
      `mu0[1]` = 0, `mu0[2]` = 0,
      `Phi0[1]` = sim$phi_fixed[1, 1], `Phi0[2]` = sim$phi_fixed[1, 2],
      `Phi0[3]` = sim$phi_fixed[2, 1], `Phi0[4]` = sim$phi_fixed[2, 2],
      rho = sim$rho
    )
    if (mod == "SNG") {
      truth <- c(truth,
        `omega[1]` = 1, `omega[2]` = 1,
        `alpha[1]` = sim$margin_info$margin1$alpha %||% 0,
        `alpha[2]` = sim$margin_info$margin2$alpha %||% 0
      )
    }

    fit <- readRDS(f)
    pars <- c(CORE, if (mod == "SNG") EXTRA)
    s <- grab(fit, pars) |>
      mutate(
        truth = truth[param],
        bias = mean - truth,
        rel_bias = ifelse(abs(truth) < .Machine$double.eps,
          bias, bias / abs(truth)
        ),
        cover95 = (l95 <= truth & u95 >= truth) * 1,
        model = mod, condition_id = cid, rep_id = rid,
        n_div = count_div(fit)
      )
  })

  write_csv(rep_tbl, file.path(res_dir, "summary_replications.csv"))

  cond_tbl <- rep_tbl |>
    group_by(condition_id, model, param) |>
    summarise(
      mean_rel_bias = mean(rel_bias),
      mean_bias = mean(bias),
      coverage_95 = mean(cover95),
      mean_post_sd = mean(sd),
      emp_sd = sd(mean),
      sd_bias = mean_post_sd - emp_sd,
      mean_n_div = mean(n_div),
      .groups = "drop"
    )
  write_csv(cond_tbl, file.path(res_dir, "summary_conditions.csv"))
  message(">>> Results saved to 'results/' folder.")
}
