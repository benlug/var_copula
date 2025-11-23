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
  real smooth_max(vector a, real kappa) {
    return log_sum_exp(kappa * a) / kappa;
  }
}
data {
  int<lower=2> T;
  matrix[T,2] y;
  int<lower=-1,upper=1> skew_direction[2];
}
parameters {
  row_vector[2] mu;
  real<lower=-0.99, upper=0.99> phi11;
  real<lower=-0.99, upper=0.99> phi12;
  real<lower=-0.99, upper=0.99> phi21;
  real<lower=-0.99, upper=0.99> phi22;
  vector[2] eta;
  real rho_raw;
}
transformed parameters {
  matrix[2,2] Phi_T;
  matrix[T-1,2] z;        // VAR residuals (innovations)
  vector[2]     sigma_exp;
  real          rho;

  rho = 0.97 * tanh(rho_raw);
  Phi_T[1,1]=phi11; Phi_T[2,1]=phi12;
  Phi_T[1,2]=phi21; Phi_T[2,2]=phi22;

  {
    matrix[T-1,2] pred = rep_matrix(mu, T-1) + y[1:T-1,] * Phi_T;
    z = y[2:T,] - pred;
  }

  {
    real kappa = 10.0;
    vector[T-1] g;
    vector[2] b;

    for (i in 1:2) {
      for (t in 1:(T-1)) g[t] = -skew_direction[i] * z[t,i];
      b[i] = smooth_max(g, kappa);              // >= max_t(-s * z)
      sigma_exp[i] = b[i] + exp(eta[i]) + 1e-9; // > 0 with epsilon
    }
  }
}
model {
  mu    ~ normal(0, 0.25);
  phi11 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  rho_raw ~ normal(0, 0.75);

  for (i in 1:2)
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];

  {
    vector[2] rate_exp;
    for (i in 1:2) rate_exp[i] = 1.0 / sigma_exp[i];

    for (t in 1:(T-1)) {
      vector[2] u;
      for (i in 1:2) {
        real x = sigma_exp[i] + skew_direction[i] * z[t,i];
        target += exponential_lpdf(x | rate_exp[i]);
        {
          real ui = exponential_cdf(x, rate_exp[i]);   // x>=0 here
          u[i] = (skew_direction[i] == 1) ? ui : (1.0 - ui);
        }
      }
      target += gaussian_copula_ld(u[1], u[2], rho);
    }
  }
}
generated quantities {
  matrix[2,2] Phi = Phi_T';
}

