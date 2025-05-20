// model_sng_sl.stan
// single-level bivariate var(1) with skew-normal margins
// dependence is modeled with a gaussian copula

functions {
  // gaussian copula log-density (same as NG model)
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9;
    real u_clamp = fmax(eps, fmin(1 - eps, u));
    real v_clamp = fmax(eps, fmin(1 - eps, v));
    if (rho <= -1.0 + eps || rho >= 1.0 - eps) return negative_infinity();
    real z1 = inv_Phi(u_clamp);
    real z2 = inv_Phi(v_clamp);
    real rho_sq = square(rho);
    if (rho_sq >= 1.0 - eps) return negative_infinity();
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

  // parameters for skew-normal margins
  vector<lower=0>[2] omega; // scale
  vector[2] alpha;          // shape

  // gaussian copula correlation
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2, 2] Phi;
  vector[2] xi;
  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;
  xi = -omega .* (alpha ./ sqrt(1 + square(alpha))) * sqrt(2 / pi());
}

model {
  // priors
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho ~ normal(0, 0.5);
  // priors for skew-normal parameters (xi determined by alpha and omega)
  omega ~ normal(0, 1);   // prior for scale (half-normal implied, ~1)
  alpha ~ cauchy(0, 5);   // reasonably wide prior for shape

  // likelihood
  for (t in 2:T) {
    vector[2] y_curr = y[t];
    vector[2] y_prev = y[t-1];
    vector[2] cond_mean = mu + Phi * y_prev;
    vector[2] residuals = y_curr - cond_mean;

    // add marginal log-likelihoods (skew-normal)
    target += skew_normal_lpdf(residuals[1] | xi[1], omega[1], alpha[1]);
    target += skew_normal_lpdf(residuals[2] | xi[2], omega[2], alpha[2]);

    // calculate cdfs and add gaussian copula density
    real u1 = skew_normal_cdf(residuals[1] | xi[1], omega[1], alpha[1]);
    real u2 = skew_normal_cdf(residuals[2] | xi[2], omega[2], alpha[2]);
    target += gaussian_copula_density(u1, u2, rho);
  }
}

generated quantities {

}
