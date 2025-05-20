// model_ng_sl.stan
// single-level bivariate var(1) with normal margins
// dependence between residuals is modeled with a gaussian copula

functions {
  // gaussian copula log-density
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9;
    real u_clamp = fmax(eps, fmin(1 - eps, u));
    real v_clamp = fmax(eps, fmin(1 - eps, v));
    if (rho <= -1.0 + eps || rho >= 1.0 - eps) return negative_infinity();
    real z1 = inv_Phi(u_clamp);
    real z2 = inv_Phi(v_clamp);
    real rho_sq = square(rho);
    if (rho_sq >= 1.0 - eps) return negative_infinity(); // avoid issues near +/- 1
    real inv_1_minus_rho_sq = 1.0 / (1.0 - rho_sq);
    return -0.5 * log1m(rho_sq)
           - 0.5 * inv_1_minus_rho_sq * (square(z1) - 2.0 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=2> T;       // number of time points
  vector[2] y[T];       // observed series y[t]
}

parameters {
  vector[2] mu;                // intercepts for y1 and y2
  real<lower=-1, upper=1> phi11; // effect of y1[t-1] on y1[t]
  real<lower=-1, upper=1> phi12; // effect of y2[t-1] on y1[t]
  real<lower=-1, upper=1> phi21; // effect of y1[t-1] on y2[t]
  real<lower=-1, upper=1> phi22; // effect of y2[t-1] on y2[t]

  // residual standard deviations for the normal margins
  vector<lower=0>[2] sigma;

  // gaussian copula correlation
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  // construct Phi matrix for convenience
  matrix[2, 2] Phi;
  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;
}

model {
  // priors (weakly informative)
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5); // centered at 0, allow moderate values
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho ~ normal(0, 0.5);   // centered at 0, allows moderate correlation
  sigma ~ normal(0, 1);   // prior on sds -> half-normal implied by lower=0 constraint

  // likelihood
  for (t in 2:T) {
    vector[2] y_curr = y[t];
    vector[2] y_prev = y[t-1];
    vector[2] cond_mean = mu + Phi * y_prev; // conditional mean
    vector[2] residuals = y_curr - cond_mean; // residuals

    // add marginal log-likelihoods (normal)
    target += normal_lpdf(residuals[1] | 0, sigma[1]);
    target += normal_lpdf(residuals[2] | 0, sigma[2]);

    // calculate cdfs and add gaussian copula density
    real u1 = normal_cdf(residuals[1] | 0, sigma[1]);
    real u2 = normal_cdf(residuals[2] | 0, sigma[2]);
    target += gaussian_copula_density(u1, u2, rho);
  }
}

generated quantities {
  // optional: calculate log-likelihood for model comparison
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
