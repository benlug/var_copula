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
  vector[2] skew_direction; // +1 for right skew, -1 for left skew
}

parameters {
  row_vector[2] mu;

  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // Unconstrained parameter for reparameterized scale
  // sigma_exp = b + exp(eta), where b is the data-dependent feasibility bound
  vector[2] eta;

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
  vector[2] sigma_exp;  // implied scale = b + exp(eta)
  vector[2] rate_exp;

  // Priors
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho   ~ normal(0, 0.5);
  
  // IMPORTANT (Study 2 fix): Do NOT place an additional independent prior on `eta`.
  // We already induce a prior on `sigma_exp` below via change-of-variables
  // (sigma_exp = b + exp(eta), Jacobian term = eta). An extra prior on eta would
  // double-regularize the scale and can strongly conflict with the intended prior
  // on sigma_exp, worsening geometry (divergences) and bias.

  // Compute residuals
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals   = y[2:T, ] - predictions;

  // Compute feasibility bounds: b_i = max_t(-s_i * res_{t,i})
  // This ensures x_shifted = sigma_exp + s_i * res > 0 for all t
  for (i in 1:2) {
    real m = -skew_direction[i] * residuals[1, i];
    for (t in 2:(T-1)) {
      m = fmax(m, -skew_direction[i] * residuals[t, i]);
    }
    b[i] = m;
  }

  // Reparameterized sigma: sigma_exp = b + exp(eta) (always feasible)
  for (i in 1:2) {
    sigma_exp[i] = b[i] + exp(eta[i]);
  }
  rate_exp = 1.0 ./ sigma_exp;

  // Induced prior on sigma_exp via change-of-variables
  // If we want sigma_exp ~ lognormal(0, 0.5), we add:
  // target += lognormal_lpdf(sigma_exp | 0, 0.5) + eta (Jacobian: d(sigma_exp)/d(eta) = exp(eta))
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
  }

  // Likelihood
  for (t in 1:(T-1)) {
    row_vector[2] res = residuals[t];
    vector[2] x_shifted;
    vector[2] u_vec;

    for (i in 1:2) {
      // Shift residual to positive support: x_shifted > 0 by construction
      x_shifted[i] = sigma_exp[i] + skew_direction[i] * res[i];
      target += exponential_lpdf(x_shifted[i] | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted[i], rate_exp[i]);
      // For left-skewed (mirrored) margins, flip the CDF
      if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
    }

    target += gaussian_copula_ld(u_vec[1], u_vec[2], rho);
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';

  // FIXED: Output names now match what analysis code expects
  vector[2] sigma_exp;  // Changed from sigma_exp_out to sigma_exp
  vector[2] b_gq;       // feasibility bounds (renamed to avoid confusion)
  vector[2] slack;      // exp(eta) = sigma_exp - b
  vector[2] rate_exp;

  {
    // Recompute residuals at the current draw
    matrix[T-1, 2] predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
    matrix[T-1, 2] residuals   = y[2:T, ] - predictions;

    for (i in 1:2) {
      real m = -skew_direction[i] * residuals[1, i];
      for (t in 2:(T-1)) m = fmax(m, -skew_direction[i] * residuals[t, i]);
      b_gq[i]      = m;
      slack[i]     = exp(eta[i]);
      sigma_exp[i] = b_gq[i] + slack[i];  // Now outputs as sigma_exp[1], sigma_exp[2]
      rate_exp[i]  = 1.0 / sigma_exp[i];
    }
  }
}
