// file: model_skewnormal_multilevel_vec.stan
// Partially vectorized version. EXPECTS y_t data format.

functions {
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
  int<lower=1> N;
  int<lower=2> T;
  // *** EXPECTS DATA AS: Array of matrices: y_t[t] is an N x 2 matrix ***
  matrix[N, 2] y_t[T];
  // Hyperparameters for priors (Only needed for multilevel)
  real<lower=0> prior_sigma_mu_scale;
  real<lower=0> prior_sigma_phi_scale;
}

parameters {
  // Hierarchical Parameters
  vector[2] mu_global;
  real<lower=-1, upper=1> phi11_base;
  real<lower=-1, upper=1> phi12_base;
  real<lower=-1, upper=1> phi21_base;
  real<lower=-1, upper=1> phi22_base;
  vector<lower=0>[2] sigma_mu;
  vector<lower=0>[4] sigma_phi;
  // Non-centered parameterization raw deviates
  vector[N] mu1_raw; vector[N] mu2_raw;
  vector[N] dev_phi11_raw; vector[N] dev_phi12_raw;
  vector[N] dev_phi21_raw; vector[N] dev_phi22_raw;

  // Residual Distribution Parameters (Skew-Normal)
  vector<lower=0>[2] omega; // scale > 0
  vector[2] alpha;          // shape
  vector[2] xi;             // location parameter

  // Copula Correlation Parameter
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  // Subject-specific parameters
  vector[N] mu1 = mu_global[1] + mu1_raw * sigma_mu[1];
  vector[N] mu2 = mu_global[2] + mu2_raw * sigma_mu[2];
  vector[N] phi11_i = phi11_base + dev_phi11_raw * sigma_phi[1];
  vector[N] phi12_i = phi12_base + dev_phi12_raw * sigma_phi[2];
  vector[N] phi21_i = phi21_base + dev_phi21_raw * sigma_phi[3];
  vector[N] phi22_i = phi22_base + dev_phi22_raw * sigma_phi[4];
}

model {
  // Priors
  mu_global ~ normal(0, 1);
  phi11_base ~ normal(0.3, 0.3); phi12_base ~ normal(0, 0.3);
  phi21_base ~ normal(0, 0.3); phi22_base ~ normal(0.3, 0.3);
  // Priors for random effect SDs
  sigma_mu ~ normal(0, prior_sigma_mu_scale);
  sigma_phi ~ normal(0, prior_sigma_phi_scale);
  // Priors for non-centered raw deviates
  mu1_raw ~ std_normal(); mu2_raw ~ std_normal();
  dev_phi11_raw ~ std_normal(); dev_phi12_raw ~ std_normal();
  dev_phi21_raw ~ std_normal(); dev_phi22_raw ~ std_normal();
  // Copula prior
  rho ~ normal(0, 0.5);
  // Priors for Skew-Normal parameters
  xi ~ normal(0, 1);      // Prior for location
  omega ~ normal(0, 1);   // Prior for scale (Half-Normal implied)
  alpha ~ cauchy(0, 5);   // Reasonably wide prior for shape, Cauchy(0,10) was quite wide

  // --- Partially Vectorized Likelihood Calculation ---
  for (t in 2:T) {
    matrix[N, 2] y_curr = y_t[t];
    matrix[N, 2] y_prev = y_t[t-1];

    // Calculate centered previous values (using subject-specific means)
    matrix[N, 2] centered_prev;
    centered_prev[, 1] = y_prev[, 1] - mu1;
    centered_prev[, 2] = y_prev[, 2] - mu2;

    // Calculate predictions (using subject-specific phis)
    vector[N] pred1 = phi11_i .* centered_prev[, 1] + phi12_i .* centered_prev[, 2];
    vector[N] pred2 = phi21_i .* centered_prev[, 1] + phi22_i .* centered_prev[, 2];

    // Calculate conditional means (using subject-specific means)
    vector[N] cond_mean1 = mu1 + pred1;
    vector[N] cond_mean2 = mu2 + pred2;

    // Calculate residuals
    vector[N] e1 = y_curr[, 1] - cond_mean1;
    vector[N] e2 = y_curr[, 2] - cond_mean2;

    // 1. Add marginal log-likelihoods (vectorized SkewNormal)
    target += skew_normal_lpdf(e1 | xi[1], omega[1], alpha[1]);
    target += skew_normal_lpdf(e2 | xi[2], omega[2], alpha[2]);

    // 2. Calculate CDFs and add Copula density (loop over subjects i)
    for (i in 1:N) {
        real u1 = skew_normal_cdf(e1[i] | xi[1], omega[1], alpha[1]);
        real u2 = skew_normal_cdf(e2[i] | xi[2], omega[2], alpha[2]);
        target += gaussian_copula_density(u1, u2, rho);
    }
  } // End loop over time t
}

generated quantities {
  real log_lik = 0; // Placeholder
  // vector[2] var_mu = square(sigma_mu);
  // vector[4] var_phi = square(sigma_phi);
}
