/*
 * model_TVP_NG_sl.stan
 * =====================
 * Time-Varying Parameter Copula-VAR(1) with:
 *   - Normal (Gaussian) marginal innovations
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
    
    // Copula log-density
    return -0.5 * log1m(rho2)
           - 0.5 / (1 - rho2) * (square(z1) - 2 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=2> T;
  matrix[T, 2] y;
}

parameters {
  // VAR parameters
  row_vector[2] mu; 
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;
  
  // Marginal innovation SDs
  vector<lower=0>[2] sigma;
  
  // State-space parameters for time-varying rho
  real z0;                           // Initial state (on Fisher-z scale)
  real<lower=0> sigma_z;             // State innovation SD
  vector[T-2] z_innov;               // Standardized state innovations (non-centered)
}

transformed parameters {
  matrix[2, 2] Phi_T;
  vector[T-1] z;                     // Latent state trajectory
  vector[T-1] rho;                   // Time-varying copula correlation
  
  // Construct Phi transpose for efficient computation
  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;
  
  // Build state trajectory (non-centered parameterization)
  z[1] = z0;
  for (t in 2:(T-1)) {
    z[t] = z[t-1] + sigma_z * z_innov[t-1];
  }
  
  // Transform to correlation scale via tanh
  for (t in 1:(T-1)) {
    rho[t] = tanh(z[t]);
  }
}

model {
  matrix[T-1, 2] residuals;
  matrix[T-1, 2] predictions;
  real log_sigma_sum = sum(log(sigma));

  // Priors on VAR parameters
  mu[1] ~ normal(0, 1);
  mu[2] ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  sigma ~ normal(0, 1);  // Half-normal (positive constraint)
  
  // Priors on state-space parameters
  z0 ~ normal(0, 1);                 // Prior on initial state
  sigma_z ~ normal(0, 0.1);          // Shrinkage prior: favors constant rho
  z_innov ~ std_normal();            // Standardized innovations
  
  // Compute residuals
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals = y[2:T, ] - predictions;

  // Likelihood with time-varying copula
  for (t in 1:(T-1)) {
    row_vector[2] res = residuals[t];
    row_vector[2] z_std = res ./ sigma';
    
    // Marginal contributions
    target += std_normal_lpdf(z_std[1]);
    target += std_normal_lpdf(z_std[2]);
    target += -log_sigma_sum;
    
    // Copula contribution with time-varying rho
    real u1 = Phi(z_std[1]);
    real u2 = Phi(z_std[2]);
    target += gaussian_copula_ld(u1, u2, rho[t]);
  }
}

generated quantities {
  matrix[2, 2] Phi = Phi_T';
  
  // Summary statistics for rho trajectory
  real rho_mean = mean(rho);
  real rho_sd = sd(rho);
  real rho_min = min(rho);
  real rho_max = max(rho);
  real rho_range = rho_max - rho_min;
  
  // First and last values
  real rho_start = rho[1];
  real rho_end = rho[T-1];
  real rho_change = rho_end - rho_start;
  
  // NOTE: Full trajectories removed to save disk space
  // If needed for diagnostics, uncomment:
  // vector[T-1] rho_trajectory = rho;
  // vector[T-1] z_trajectory = z;
}
