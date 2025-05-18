// file: model_skewnormal_singlelevel_vec.stan
// Non-hierarchical (single-level) vectorized SkewNormal margins model.
// v2: Removed incorrect is_finite check. Added clamping before copula call.

functions {
  // Gaussian copula log-density
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9; // Small epsilon
    // Clamp inputs to copula density to be safe
    real u_clamp = fmax(eps, fmin(1.0 - eps, u));
    real v_clamp = fmax(eps, fmin(1.0 - eps, v));

    if (rho <= -1.0 + eps || rho >= 1.0 - eps) return negative_infinity();

    real z1 = inv_Phi(u_clamp);
    real z2 = inv_Phi(v_clamp);
    // Removed the incorrect is_finite check here

    real rho_sq = square(rho);
    // Check rho_sq to avoid issues in log1m and division
    if (rho_sq >= 1.0 - eps) return negative_infinity(); // Check if rho is too close to +/- 1

    real inv_1_minus_rho_sq = 1.0 / (1.0 - rho_sq);
    return -0.5 * log1m(rho_sq) // Use log1m for potentially better precision
           - 0.5 * inv_1_minus_rho_sq * (square(z1) - 2.0 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=1> N;
  int<lower=2> T;
  matrix[N, 2] y_t[T];
  // NO Prior scales for RE needed
}

parameters {
  // Global Parameters
  vector[2] mu_global;
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // Residual Distribution Parameters (Skew-Normal)
  vector<lower=0>[2] omega; // scale > 0
  vector[2] alpha;          // shape
  vector[2] xi;             // location

  // Copula Correlation Parameter
  real<lower=-1, upper=1> rho;
}

// NO transformed parameters block for random effects

model {
  // Priors
  mu_global ~ normal(0, 1);
  phi11 ~ normal(0.3, 0.3); phi12 ~ normal(0, 0.3);
  phi21 ~ normal(0, 0.3); phi22 ~ normal(0.3, 0.3);
  rho ~ normal(0, 0.5);
  xi ~ normal(0, 1);
  omega ~ normal(0, 1); // Half-Normal implied by lower=0 constraint
  alpha ~ cauchy(0, 5);

  // --- Vectorized Likelihood Calculation ---
  for (t in 2:T) {
    matrix[N, 2] y_curr = y_t[t];
    matrix[N, 2] y_prev = y_t[t-1];

    matrix[N, 2] centered_prev;
    centered_prev[, 1] = y_prev[, 1] - mu_global[1];
    centered_prev[, 2] = y_prev[, 2] - mu_global[2];

    vector[N] pred1 = phi11 * centered_prev[, 1] + phi12 * centered_prev[, 2];
    vector[N] pred2 = phi21 * centered_prev[, 1] + phi22 * centered_prev[, 2];

    vector[N] cond_mean1 = mu_global[1] + pred1;
    vector[N] cond_mean2 = mu_global[2] + pred2;

    vector[N] e1 = y_curr[, 1] - cond_mean1;
    vector[N] e2 = y_curr[, 2] - cond_mean2;

    // 1. Add marginal log-likelihoods (vectorized SkewNormal)
    target += skew_normal_lpdf(e1 | xi[1], omega[1], alpha[1]);
    target += skew_normal_lpdf(e2 | xi[2], omega[2], alpha[2]);

    // 2. Calculate CDFs and add Copula density (loop over subjects i)
    real eps_cdf = 1e-9; // Epsilon for CDF clamping
    for (i in 1:N) {
        real u1_raw = skew_normal_cdf(e1[i] | xi[1], omega[1], alpha[1]);
        real u2_raw = skew_normal_cdf(e2[i] | xi[2], omega[2], alpha[2]);

        // *** Clamp CDF results HERE before passing to copula ***
        real u1 = fmax(eps_cdf, fmin(1.0 - eps_cdf, u1_raw));
        real u2 = fmax(eps_cdf, fmin(1.0 - eps_cdf, u2_raw));

        // Add copula density contribution
        target += gaussian_copula_density(u1, u2, rho);
    }
  } // End loop over time t
}

generated quantities {
  real log_lik = 0; // Placeholder
  // Note: No random effect variances to calculate here
}
