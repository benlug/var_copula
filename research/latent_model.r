library(MASS) # for mvrnorm

# ---- Utility to sample from a bivariate Clayton copula with normal(0, sd_j) marginals
rClaytonNormal <- function(n, alpha, sd = c(1, 1)) {
  # 1) sample (U1,U2) from bivariate Clayton(alpha>0)
  #    F(U1,U2) = (U1^(-alpha)+U2^(-alpha)-1)^(-1/alpha)
  #    We'll do a standard approach:
  u1 <- runif(n)
  v <- runif(n)
  # conditional formula:
  #   U2|U1 = (v*(u1^(-alpha)-1)+1)^(-1/alpha)
  u2 <- (v * (u1^(-alpha) - 1) + 1)^(-1 / alpha)

  # 2) map U1,U2 -> normal(0, sd_j)
  Z1 <- qnorm(u1, 0, sd[1])
  Z2 <- qnorm(u2, 0, sd[2])

  cbind(Z1, Z2)
}


# ---- Simulation of latent VAR(1) with either Gaussian or Clayton copula for residuals
sim_biv_latentVAR <- function(
    T = 200,
    phi = matrix(c(0.6, 0.2, 0.1, 0.5), 2, 2, byrow = TRUE),
    mu = c(0, 0),
    # copula choice & params
    copula = c("gaussian", "clayton"),
    corr_gauss = 0.5, # correlation if gaussian
    clayton_alpha = 1, # if clayton
    sigmas = c(1, 1), # marginal sds
    # thresholds for var2 -> ordinal
    thresholds = c(-0.5, 0.5)) {
  copula <- match.arg(copula)
  Z <- matrix(0, nrow = T, ncol = 2)
  # sample residuals e_t:
  if (copula == "gaussian") {
    # bivariate normal with correlation corr_gauss
    Sig <- matrix(c(
      sigmas[1]^2, corr_gauss * prod(sigmas),
      corr_gauss * prod(sigmas), sigmas[2]^2
    ), 2, 2)
    E <- mvrnorm(T, mu = c(0, 0), Sigma = Sig)
  } else {
    # clayton
    E <- rClaytonNormal(T, alpha = clayton_alpha, sd = sigmas)
  }

  # generate Z
  Z[1, ] <- mu + E[1, ] # initial
  for (t in 2:T) {
    Z[t, ] <- mu + phi %*% (Z[t - 1, ] - mu) + E[t, ]
  }

  # Observed Y1 = Z1 (continuous)
  Y1 <- Z[, 1]
  # Observed Y2 = ordinal
  # cat1 if Z2 < thresholds[1], cat2 if in [t1, t2), cat3 if >= t2
  Y2 <- cut(Z[, 2],
    breaks = c(-Inf, thresholds, Inf),
    labels = FALSE
  ) # 1,2,3

  data.frame(t = 1:T, Z1 = Z[, 1], Z2 = Z[, 2], Y1 = Y1, Y2 = Y2)
}

# example usage
set.seed(101)
dat_gauss <- sim_biv_latentVAR(T = 300, copula = "gaussian")
dat_clay <- sim_biv_latentVAR(T = 300, copula = "clayton", clayton_alpha = 1.5)




library(this.path)
setwd(this.dir())
library(rstan)
rstan_options(auto_write = TRUE)

# 1. Simulate some data from the GAUSSIAN copula for demonstration
set.seed(123)
dat_sim <- sim_biv_latentVAR(
  T = 200, copula = "gaussian",
  corr_gauss = 0.4
) # true correlation=0.4

# Suppose the thresholds for Y2 are c(-0.5, 0.5) => 3 categories
K <- 4
thresh <- c(-1, 0, 1)

stan_data <- list(
  T = nrow(dat_sim),
  Y1 = dat_sim$Y1,
  Y2 = dat_sim$Y2,
  K = K,
  thresh = thresh
)


stan_code <- "
  data {
    int<lower=2> T;              // number of time points
    vector[T] Y1;                // observed continuous
    int<lower=1> Y2[T];          // observed ordinal, e.g. 1..K
    int<lower=2> K;              // number of categories for Y2
    vector[K-1] thresh;          // known thresholds for ordinal
  }
  parameters {
    vector[2] mu;                // latent means
    matrix[2,2] Phi;             // VAR(1) coefficients
    cholesky_factor_cov[2] Lres; // Cholesky factor for residual cov (2x2)
    real<lower=0,upper=1e-1> sigma_meas; // small measurement error for Y1

    // latent states Z[t] = (Z1[t], Z2[t])
    // We separate them out for clarity
    vector[T] z1;
    vector[T] z2;
  }
  transformed parameters {
    cov_matrix[2] Sigma;
    Sigma = multiply_lower_tri_self_transpose(Lres); // Lres*Lres^T
  }
  model {
    // Priors (very basic)
    mu ~ normal(0, 2);
    to_vector(Phi) ~ normal(0, 0.5);  // fairly loose
    // enforce stationarity? (omitted for brevity, or we do a prior on spectral radius)
    Lres ~ lkj_corr_cholesky(2);

    // small meas error
    sigma_meas ~ normal(0, 0.1);

    // Likelihood:
    // 1) define the dynamic for Z[t], t=2..T
    //    Z[t] ~ Normal( mu + Phi*(Z[t-1] - mu), Sigma )
    for(t in 2:T) {
      vector[2] mean_t;
      mean_t[1] = mu[1] + Phi[1,1]*(z1[t-1]-mu[1]) + Phi[1,2]*(z2[t-1]-mu[2]);
      mean_t[2] = mu[2] + Phi[2,1]*(z1[t-1]-mu[1]) + Phi[2,2]*(z2[t-1]-mu[2]);

      // multi_normal: we store (z1[t], z2[t]) in a vector
      vector[2] zt;
      zt[1] = z1[t];
      zt[2] = z2[t];
      zt ~ multi_normal_cholesky(mean_t, Lres);
    }

    // 2) tie (z1[t]) to observed Y1[t] with small measurement error
    //    => Y1[t] ~ Normal(z1[t], sigma_meas)
    for(t in 1:T) {
      Y1[t] ~ normal(z1[t], sigma_meas);
    }

    // 3) ordinal measurement for Y2:
    //    if Y2[t] == 1 => z2[t] < thresh[1]
    //    if Y2[t] == 2 => thresh[1] <= z2[t] < thresh[2]
    //    ...
    //    if Y2[t] == K => z2[t] >= thresh[K-1]
    for(t in 1:T) {
      if(Y2[t] == 1) {
        target += normal_lcdf(thresh[1] | z2[t], 1e6) // hack: no scale
                  - normal_lcdf(-1e9 | z2[t], 1e6);   // effectively z2 < thresh[1]
      } else if(Y2[t] == K) {
        target += normal_lcdf(1e9 | z2[t], 1e6)
                  - normal_lcdf(thresh[K-1] | z2[t], 1e6);
      } else {
        // category c => thresh[c-1] <= z2[t] < thresh[c]
        int c = Y2[t];
        target += normal_lcdf(thresh[c] | z2[t], 1e6)
                  - normal_lcdf(thresh[c-1] | z2[t], 1e6);
      }
    }
  }
"

# 2. Compile and run
mod <- stan_model(model_code = stan_code)
fit <- sampling(mod, data = stan_data, chains = 4, iter = 4000, warmup = 2000, cores = parallel::detectCores(), control = list(adapt_delta = 0.95, max_treedepth = 12))

print(fit, pars = c("mu", "Phi", "Sigma", "sigma_meas"), probs = c(0.05, 0.95))
