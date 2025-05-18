// file: model_normal_singlelevel_vec.stan
// Non-hierarchical (single-level) vectorized Normal margins model.
// Assumes all subjects share the same intercept and VAR parameters.

functions {
  // Gaussian copula log-density (same as before)
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9;
    real u_clamp = fmax(eps, fmin(1 - eps, u));
    real v_clamp = fmax(eps, fmin(1 - eps, v));
    if (rho <= -1.0 + eps || rho >= 1.0 - eps) return negative_infinity();
    real z1 = inv_Phi(u_clamp);
    real z2 = inv_Phi(v_clamp);
    real rho_sq = square(rho);
    if (rho_sq >= 1.0) return negative_infinity();
    real inv_1_minus_rho_sq = 1.0 / (1.0 - rho_sq);
    return -0.5 * log1m(rho_sq)
           - 0.5 * inv_1_minus_rho_sq * (square(z1) - 2.0 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=1> N; // Number of subjects (still needed for data shape)
  int<lower=2> T;
  matrix[N, 2] y_t[T]; // Data remains N x 2 per time point
  // NO Prior scales for RE needed in data block
}

parameters {
  // Global Parameters Only
  vector[2] mu_global;          // Global intercept
  real<lower=-1, upper=1> phi11; // Global VAR param
  real<lower=-1, upper=1> phi12; // Global VAR param
  real<lower=-1, upper=1> phi21; // Global VAR param
  real<lower=-1, upper=1> phi22; // Global VAR param

  // Residual Distribution Parameters (Normal)
  vector<lower=0>[2] sigma;      // Residual SDs

  // Copula Correlation Parameter
  real<lower=-1, upper=1> rho;
}

// NO transformed parameters block for random effects

model {
  // Priors (only on global parameters)
  mu_global ~ normal(0, 1);
  phi11 ~ normal(0.3, 0.3); phi12 ~ normal(0, 0.3);
  phi21 ~ normal(0, 0.3); phi22 ~ normal(0.3, 0.3);
  rho ~ normal(0, 0.5);
  sigma ~ normal(0, 1); // Prior on [sigma1, sigma2] -> Half-Normal implied

  // --- Vectorized Likelihood Calculation (Simpler) ---
  for (t in 2:T) {
    matrix[N, 2] y_curr = y_t[t];
    matrix[N, 2] y_prev = y_t[t-1];

    // Calculate centered previous values (vectorized, using global mean)
    matrix[N, 2] centered_prev;
    centered_prev[, 1] = y_prev[, 1] - mu_global[1];
    centered_prev[, 2] = y_prev[, 2] - mu_global[2];

    // Calculate predictions (vectorized, using global phis)
    vector[N] pred1 = phi11 * centered_prev[, 1] + phi12 * centered_prev[, 2];
    vector[N] pred2 = phi21 * centered_prev[, 1] + phi22 * centered_prev[, 2];

    // Calculate conditional means (vectorized)
    vector[N] cond_mean1 = mu_global[1] + pred1;
    vector[N] cond_mean2 = mu_global[2] + pred2;

    // Calculate residuals (vectorized)
    vector[N] e1 = y_curr[, 1] - cond_mean1;
    vector[N] e2 = y_curr[, 2] - cond_mean2;

    // 1. Add marginal log-likelihoods (vectorized Normal)
    target += normal_lpdf(e1 | 0, sigma[1]);
    target += normal_lpdf(e2 | 0, sigma[2]);

    // 2. Calculate CDFs and add Copula density (loop over subjects i)
    for (i in 1:N) {
        real u1 = normal_cdf(e1[i] | 0, sigma[1]);
        real u2 = normal_cdf(e2[i] | 0, sigma[2]);
        target += gaussian_copula_density(u1, u2, rho);
    }
  } // End loop over time t
}

generated quantities {
  real log_lik = 0; // Placeholder
  // Note: No random effect variances to calculate here
}
