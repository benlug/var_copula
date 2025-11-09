// Study A: indicator-skew (measurement residuals), Gaussian state.
// Non-centered state path; feasibility-lift uses e_t = y_t - state_t.

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
  int<lower=1> T;
  matrix[T,2] y;
  int<lower=-1,upper=1> skew_direction[2];
}

parameters {
  vector[2] mu;

  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  vector[2] eta;                      // sigma_exp = b + exp(eta)
  real<lower=-0.995, upper=0.995> rho;

  // Non-centered state innovations
  matrix[T,2] z_raw;                  // ~ N(0,1)
}

transformed parameters {
  matrix[2,2] B;
  matrix[T,2] state;

  B[1,1]=phi11; B[1,2]=phi12;
  B[2,1]=phi21; B[2,2]=phi22;

  {
    vector[2] s;
    for (t in 1:T) {
      vector[2] zt = to_vector(z_raw[t]);
      if (t == 1) s = mu + zt;
      else        s = mu + B * s + zt;
      state[t,1] = s[1];
      state[t,2] = s[2];
    }
  }
}

model {
  // ----- Priors -----
  mu    ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.3);
  phi21 ~ normal(0, 0.3);
  to_vector(z_raw) ~ normal(0, 1);
  rho   ~ normal(0, 0.5);

  // ----- Measurement residuals and feasibility bounds -----
  matrix[T,2] e;
  for (t in 1:T) {
    e[t,1] = y[t,1] - state[t,1];
    e[t,2] = y[t,2] - state[t,2];
  }

  vector[2] b; // b_i = max_t{-d_i * e_{t,i}}
  for (i in 1:2) {
    real m = -skew_direction[i] * e[1,i];
    for (t in 2:T) m = fmax(m, -skew_direction[i] * e[t,i]);
    b[i] = m;
  }

  vector[2] sigma_exp = b + exp(eta);
  vector[2] rate_exp  = 1.0 ./ sigma_exp;

  // Prior on sigma_exp via change-of-variables
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
  }

  // ----- Likelihood -----
  for (t in 1:T) {
    vector[2] x;
    vector[2] u;
    for (i in 1:2) {
      x[i] = sigma_exp[i] + skew_direction[i] * e[t,i]; // > 0 by construction
      target += exponential_lpdf(x[i] | rate_exp[i]);
      u[i]    = exponential_cdf(x[i], rate_exp[i]);
      if (skew_direction[i] < 0) u[i] = 1.0 - u[i];
    }
    target += gaussian_copula_ld(u[1], u[2], rho);
  }
}

generated quantities {
  // Export B if you want it for diagnostics
}
