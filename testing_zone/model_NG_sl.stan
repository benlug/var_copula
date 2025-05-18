// file: model_NG_sl.stan
// Bivariate VAR(1) model assuming Normal margins and Gaussian copula.
// Single-level (no random effects).

functions {
  // Gaussian copula log-density
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9;
    real u_clamp = fmax(eps, fmin(1 - eps, u));
    real v_clamp = fmax(eps, fmin(1 - eps, v));
    if (rho <= -1.0 + eps || rho >= 1.0 - eps) return negative_infinity();
    real z1 = inv_Phi(u_clamp);
    real z2 = inv_Phi(v_clamp);
    real rho_sq = square(rho);
    if (rho_sq >= 1.0 - eps) return negative_infinity(); // Avoid issues near +/- 1
    real inv_1_minus_rho_sq = 1.0 / (1.0 - rho_sq);
    return -0.5 * log1m(rho_sq)
           - 0.5 * inv_1_minus_rho_sq * (square(z1) - 2.0 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=2> T;       // Number of time points
  vector[2] y[T];       // Data: array of 2-element vectors y_t
}

parameters {
  // Global Parameters Only
  vector[2] mu;                // Intercepts (fixed at 0 in sim, but estimate)
  real<lower=-1, upper=1> phi11; // VAR param
  real<lower=-1, upper=1> phi12; // VAR param
  real<lower=-1, upper=1> phi21; // VAR param
  real<lower=-1, upper=1> phi22; // VAR param

  // Residual Distribution Parameters (Normal)
  vector<lower=0>[2] sigma;      // Residual SDs (sigma1, sigma2)

  // Copula Correlation Parameter
  real<lower=-1, upper=1> rho;   // Gaussian copula correlation
}

transformed parameters {
  // Construct Phi matrix for convenience (optional)
  matrix[2, 2] Phi;
  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;
}

model {
  // Priors (Weakly informative)
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5); // Centered at 0, allow moderate values
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho ~ normal(0, 0.5);   // Centered at 0, allows moderate correlation
  sigma ~ normal(0, 1);   // Prior on SDs -> Half-Normal implied by lower=0 constraint

  // Likelihood Calculation
  for (t in 2:T) {
    vector[2] y_curr = y[t];
    vector[2] y_prev = y[t-1];
    vector[2] cond_mean = mu + Phi * y_prev; // Calculate conditional mean
    vector[2] residuals = y_curr - cond_mean; // Calculate residuals

    // 1. Add marginal log-likelihoods (Normal)
    target += normal_lpdf(residuals[1] | 0, sigma[1]);
    target += normal_lpdf(residuals[2] | 0, sigma[2]);

    // 2. Calculate CDFs and add Gaussian Copula density
    real u1 = normal_cdf(residuals[1] | 0, sigma[1]);
    real u2 = normal_cdf(residuals[2] | 0, sigma[2]);
    target += gaussian_copula_density(u1, u2, rho);
  }
}

generated quantities {
  // Optional: Calculate log-likelihood if needed for model comparison
  // real log_lik = 0;
  // for (t in 2:T) {
  //   vector[2] y_curr = y[t];
  //   vector[2] y_prev = y[t-1];
  //   vector[2] cond_mean = mu + Phi * y_prev;
  //   vector[2] residuals = y_curr - cond_mean;
  //   log_lik += normal_lpdf(residuals[1] | 0, sigma[1]);
  //   log_lik += normal_lpdf(residuals[2] | 0, sigma[2]);
  //   real u1 = normal_cdf(residuals[1] | 0, sigma[1]);
  //   real u2 = normal_cdf(residuals[2] | 0, sigma[2]);
  //   log_lik += gaussian_copula_density(u1, u2, rho);
  // }
}
