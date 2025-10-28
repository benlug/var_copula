functions {
  real gaussian_copula_ld(real u, real v, real rho) {
    real eps = 1e-9;
    real uu  = fmax(eps, fmin(1 - eps, u));
    real vv  = fmax(eps, fmin(1 - eps, v));
    real z1  = inv_Phi(uu);
    real z2  = inv_Phi(vv);
    real rho2 = square(rho);
    return -0.5 * log1m(rho2)
           - 0.5 / (1 - rho2) * (square(z1) - 2 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=2> T;
  matrix[T, 2] y;
  vector[2] skew_direction; // +/-1
}

parameters {
  row_vector[2] mu;

  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // NEW: unconstrained excess over the feasibility bound
  vector[2] eta;  // sigma_exp = b + exp(eta)

  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  Phi_T[1,1] = phi11;  Phi_T[2,1] = phi12;
  Phi_T[1,2] = phi21;  Phi_T[2,2] = phi22;
}

model {
  matrix[T-1, 2] predictions;
  matrix[T-1, 2] residuals;
  vector[2] b;          // feasibility bounds b_i(mu, Phi, y)
  vector[2] sigma_exp;  // implied scale/shift
  vector[2] rate_exp;

  // Priors (as before)
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho   ~ normal(0, 0.5);

  // Residuals
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals   = y[2:T, ] - predictions;

  // Compute feasibility bounds b_i = max_t(- s_i * res_{t,i})
  for (i in 1:2) {
    real m = -skew_direction[i] * residuals[1, i];
    for (t in 2:(T-1)) {
      m = fmax(m, -skew_direction[i] * residuals[t, i]);
    }
    b[i] = m; // may be near 1 with standardized exp
  }

  // Reparameterized sigma: sigma = b + exp(eta)  (always feasible)
  for (i in 1:2) {
    sigma_exp[i] = b[i] + exp(eta[i]);
  }
  rate_exp = 1.0 ./ sigma_exp;

  // Prior on sigma_exp with change-of-variables adjustment
  // (keeps sigma_exp ~ lognormal(0, 0.5))
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
  }

  // Likelihood
  for (t in 1:(T-1)) {
    row_vector[2] res = residuals[t];
    vector[2] x_shifted;
    vector[2] u_vec;

    for (i in 1:2) {
      x_shifted[i] = sigma_exp[i] + skew_direction[i] * res[i]; // > 0 by construction
      target += exponential_lpdf(x_shifted[i] | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted[i], rate_exp[i]);
      if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
    }

    target += gaussian_copula_ld(u_vec[1], u_vec[2], rho);
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';

  vector[2] b_out;            // feasibility bounds b_i(mu, Phi, y)
  vector[2] slack_out;        // exp(eta) = sigma - b
  vector[2] sigma_exp_out;    // implied scale = b + exp(eta)
  vector[2] rate_exp_out;     // implied rate  = 1 / sigma

  {
    // Recompute residuals at the current draw
    matrix[T-1, 2] predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
    matrix[T-1, 2] residuals   = y[2:T, ] - predictions;

    for (i in 1:2) {
      real m = -skew_direction[i] * residuals[1, i];
      for (t in 2:(T-1)) m = fmax(m, -skew_direction[i] * residuals[t, i]);
      b_out[i]         = m;
      slack_out[i]     = exp(eta[i]);
      sigma_exp_out[i] = b_out[i] + slack_out[i];
      rate_exp_out[i]  = 1.0 / sigma_exp_out[i];
    }
  }
}
