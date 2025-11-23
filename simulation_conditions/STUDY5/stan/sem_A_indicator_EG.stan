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
  // Smooth upper bound of max(a): (1/kappa) * log_sum_exp(kappa * a), kappa > 0
  real smooth_max(vector a, real kappa) {
    return log_sum_exp(kappa * a) / kappa;
  }
}
data {
  int<lower=1> T;
  matrix[T,2] y;
  int<lower=-1,upper=1> skew_direction[2]; // s = +/-1
}
parameters {
  vector[2] mu;
  real<lower=-0.99, upper=0.99> phi11;
  real<lower=-0.99, upper=0.99> phi12;
  real<lower=-0.99, upper=0.99> phi21;
  real<lower=-0.99, upper=0.99> phi22;
  vector[2] eta;       // log-slack in sigma = b + exp(eta)
  real rho_raw;        // unconstrained; rho = 0.97 * tanh(rho_raw)
  matrix[T,2] z_raw;   // Gaussian state innovations
}
transformed parameters {
  matrix[2,2] B;
  matrix[T,2] state;
  matrix[T,2] e;
  vector[2]   sigma_exp;
  real        rho;

  rho = 0.97 * tanh(rho_raw); // keeps |rho| < 0.97 away from singular edge
  B[1,1]=phi11; B[1,2]=phi12; B[2,1]=phi21; B[2,2]=phi22;

  { // latent state recursion (non-centered)
    vector[2] s;
    for (t in 1:T) {
      vector[2] zt = to_vector(z_raw[t]);
      s = (t==1) ? (mu + zt) : (mu + B * s + zt);
      state[t,1] = s[1];
      state[t,2] = s[2];
    }
  }

  // measurement residuals
  for (t in 1:T) {
    e[t,1] = y[t,1] - state[t,1];
    e[t,2] = y[t,2] - state[t,2];
  }

  // smooth feasibility-lifted scales (kappa tunes sharpness)
  {
    real kappa = 10.0;
    vector[T] g;
    vector[2] b;

    for (i in 1:2) {
      for (t in 1:T) g[t] = -skew_direction[i] * e[t,i];
      b[i] = smooth_max(g, kappa);              // >= max_t(-s * e)
      sigma_exp[i] = b[i] + exp(eta[i]) + 1e-9; // > 0 with epsilon
    }
  }
}
model {
  // mild centering priors (match the Study brief)
  mu    ~ normal(0, 0.25);
  phi11 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.3);
  phi21 ~ normal(0, 0.3);
  rho_raw ~ normal(0, 0.75);
  to_vector(z_raw) ~ normal(0, 1);

  // prior on sigma with change-of-variables (sigma = b + exp(eta))
  for (i in 1:2)
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];

  // measurement-layer likelihood with Exp margins + Gaussian copula
  {
    vector[2] rate_exp;
    for (i in 1:2) rate_exp[i] = 1.0 / sigma_exp[i];

    for (t in 1:T) {
      vector[2] u;
      for (i in 1:2) {
        real x = sigma_exp[i] + skew_direction[i] * e[t,i]; // > 0 by construction
        target += exponential_lpdf(x | rate_exp[i]);
        {
          real ui = exponential_cdf(x, rate_exp[i]);        // x>=0 here
          u[i] = (skew_direction[i] == 1) ? ui : (1.0 - ui);
        }
      }
      target += gaussian_copula_ld(u[1], u[2], rho);
    }
  }
}
generated quantities {
  matrix[2,2] Phi;
  Phi[1,1]=phi11; Phi[1,2]=phi12; Phi[2,1]=phi21; Phi[2,2]=phi22;
}

