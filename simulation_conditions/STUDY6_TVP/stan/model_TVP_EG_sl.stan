/*
 * model_TVP_EG_sl.stan
 * ====================
 * Time-Varying Parameter Copula-VAR(1) with:
 *   - Exponential marginal innovations (with boundary-aware reparameterization)
 *   - Gaussian copula with time-varying correlation rho_t
 *   - State-space dynamics: z_t = z_{t-1} + eta_t, rho_t = tanh(z_t)
 *
 * Study 6: Detecting Dynamic Dependence in Intensive Longitudinal Data
 */

functions {
  /**
   * Log density of the bivariate Gaussian copula.
   */
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
  vector[2] skew_direction;  // +1 for right skew, -1 for left skew
}

parameters {
  // VAR parameters
  row_vector[2] mu;
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;
  
  // Unconstrained parameter for exponential scale reparameterization
  vector[2] eta;
  
  // State-space parameters for time-varying rho
  real z0;                           // Initial state (Fisher-z scale)
  real<lower=0> sigma_z;             // State innovation SD
  vector[T-2] z_innov;               // Standardized state innovations
}

transformed parameters {
  matrix[2, 2] Phi_T;
  vector[T-1] z;                     // Latent state trajectory
  vector[T-1] rho;                   // Time-varying copula correlation
  
  // Construct Phi transpose
  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;
  
  // Build state trajectory (non-centered parameterization)
  z[1] = z0;
  for (t in 2:(T-1)) {
    z[t] = z[t-1] + sigma_z * z_innov[t-1];
  }
  
  // Transform to correlation scale
  for (t in 1:(T-1)) {
    rho[t] = tanh(z[t]);
  }
}

model {
  matrix[T-1, 2] predictions;
  matrix[T-1, 2] residuals;
  vector[2] b;          // Feasibility bounds
  vector[2] sigma_exp;  // Implied scale = b + exp(eta)
  vector[2] rate_exp;

  // Priors on VAR parameters
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  
  // Prior on eta
  eta ~ normal(0, 1);
  
  // Priors on state-space parameters
  z0 ~ normal(0, 1);
  sigma_z ~ normal(0, 0.1);          // Shrinkage prior
  z_innov ~ std_normal();

  // Compute residuals
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals = y[2:T, ] - predictions;

  // Compute feasibility bounds
  for (i in 1:2) {
    real m = -skew_direction[i] * residuals[1, i];
    for (t in 2:(T-1)) {
      m = fmax(m, -skew_direction[i] * residuals[t, i]);
    }
    b[i] = m;
  }

  // Reparameterized sigma
  for (i in 1:2) {
    sigma_exp[i] = b[i] + exp(eta[i]);
  }
  rate_exp = 1.0 ./ sigma_exp;

  // Induced prior on sigma_exp
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
  }

  // Likelihood with time-varying copula
  for (t in 1:(T-1)) {
    row_vector[2] res = residuals[t];
    vector[2] x_shifted;
    vector[2] u_vec;

    for (i in 1:2) {
      x_shifted[i] = sigma_exp[i] + skew_direction[i] * res[i];
      target += exponential_lpdf(x_shifted[i] | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted[i], rate_exp[i]);
      if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
    }

    target += gaussian_copula_ld(u_vec[1], u_vec[2], rho[t]);
  }
}

generated quantities {
  matrix[2, 2] Phi = Phi_T';
  
  // Exponential scale parameters
  vector[2] sigma_exp;
  vector[2] b_gq;
  vector[2] slack;
  vector[2] rate_exp;

  {
    matrix[T-1, 2] predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
    matrix[T-1, 2] residuals = y[2:T, ] - predictions;

    for (i in 1:2) {
      real m = -skew_direction[i] * residuals[1, i];
      for (t in 2:(T-1)) m = fmax(m, -skew_direction[i] * residuals[t, i]);
      b_gq[i] = m;
      slack[i] = exp(eta[i]);
      sigma_exp[i] = b_gq[i] + slack[i];
      rate_exp[i] = 1.0 / sigma_exp[i];
    }
  }
  
  // Summary statistics for rho trajectory
  real rho_mean = mean(rho);
  real rho_sd = sd(rho);
  real rho_min = min(rho);
  real rho_max = max(rho);
  real rho_range = rho_max - rho_min;
  
  real rho_start = rho[1];
  real rho_end = rho[T-1];
  real rho_change = rho_end - rho_start;
  
  // NOTE: Full trajectories removed to save disk space
  // If needed for diagnostics, uncomment:
  // vector[T-1] rho_trajectory = rho;
  // vector[T-1] z_trajectory = z;
}
