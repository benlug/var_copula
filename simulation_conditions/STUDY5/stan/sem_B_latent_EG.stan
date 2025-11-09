// Study B: latent-skew (innovations), no measurement error.
// Feasibility-lift: sigma_exp = b + exp(eta), with
// b_i = max_t{-d_i * z_{t,i}} built from residuals z_t.

functions {
  real gaussian_copula_ld(real u, real v, real rho) {
    real eps = 1e-9;
    real uu  = fmax(eps, fmin(1 - eps, u));
    real vv  = fmax(eps, fmin(1 - eps, v));
    real z1  = inv_Phi(uu);
    real z2  = inv_Phi(vv);
    real r2  = square(rho);
    return -0.5 * log1m(r2)
         - 0.5 / (1 - r2) * (square(z1) - 2 * rho * z1 * z2 + square(z2))
         + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=2> T;
  matrix[T,2] y;
  int<lower=-1,upper=1> skew_direction[2]; // +1 right, -1 left
}

parameters {
  row_vector[2] mu;

  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  vector[2] eta; // sigma_exp = b + exp(eta)

  real<lower=-0.995, upper=0.995> rho;
}

transformed parameters {
  // Note: we store Phi^T so that y * Phi_T is row-wise VAR prediction
  matrix[2,2] Phi_T;
  Phi_T[1,1] = phi11;  Phi_T[2,1] = phi12;
  Phi_T[1,2] = phi21;  Phi_T[2,2] = phi22;
}

model {
  // ----- Priors -----
  mu    ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.3);
  phi21 ~ normal(0, 0.3);
  rho   ~ normal(0, 0.5);

  // ----- Residuals z_t and feasibility bounds b_i -----
  matrix[T-1,2] pred = rep_matrix(mu, T-1) + y[1:T-1,] * Phi_T;
  matrix[T-1,2] z    = y[2:T,] - pred;

  vector[2] b; // b_i = max_t{-d_i * z_{t,i}}
  for (i in 1:2) {
    real m = -skew_direction[i] * z[1,i];
    for (t in 2:(T-1)) m = fmax(m, -skew_direction[i] * z[t,i]);
    b[i] = m;
  }

  // Feasible sigma_exp and implied rates
  vector[2] sigma_exp = b + exp(eta);
  vector[2] rate_exp  = 1.0 ./ sigma_exp;

  // Prior on sigma_exp via CoV adjustment (keeps sigma roughly lognormal)
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
  }

  // ----- Likelihood with shifted exponentials + Gaussian copula -----
  for (t in 1:(T-1)) {
    vector[2] x;   // shifted argument > 0
    vector[2] u;   // PIT uniforms

    for (i in 1:2) {
      x[i] = sigma_exp[i] + skew_direction[i] * z[t,i]; // strictly > 0
      target += exponential_lpdf(x[i] | rate_exp[i]);
      u[i] = exponential_cdf(x[i], rate_exp[i]);
      if (skew_direction[i] < 0) u[i] = 1.0 - u[i];
    }
    target += gaussian_copula_ld(u[1], u[2], rho);
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';

  // Optional: export b, slack, sigma for diagnostics
  // (recompute as in model block)
  // Not shown to keep GQ short.
}
