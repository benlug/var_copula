functions {
  /**
   * Log density of the bivariate Gaussian copula.
   * Renamed from _lpdf to _ld to avoid Stan's syntax requirements for distributions.
   */
  real gaussian_copula_ld(real u, real v, real rho) {
    real eps = 1e-9;
    // Clamping prevents inv_Phi returning infinity at 0 or 1
    real uu  = fmax(eps, fmin(1 - eps, u));
    real vv  = fmax(eps, fmin(1 - eps, v));
    
    real z1  = inv_Phi(uu);
    real z2  = inv_Phi(vv);
    real rho2 = square(rho);
    
    // Log density calculation: log(c(u,v))
    // Derived from MVN PDF minus marginal PDFs
    return -0.5 * log1m(rho2)
           - 0.5 / (1 - rho2) * (square(z1) - 2 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=2> T;
  // Input format updated to matrix for efficiency
  matrix[T, 2] y;
}

parameters {
  // Defined as row_vector for efficient broadcasting in vectorized operations
  row_vector[2] mu; 
  
  // AR parameters
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;
  
  // Innovation SDs and correlation
  vector<lower=0>[2] sigma;
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  // Transposed Phi matrix for efficient matrix multiplication: Y_lag * Phi'
  matrix[2, 2] Phi_T;
  
  // Construct the transpose
  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;
}

model {
  matrix[T-1, 2] residuals;
  matrix[T-1, 2] predictions;
  
  // Pre-calculate Jacobian adjustment for standardization: log(sigma1) + log(sigma2)
  real log_sigma_sum = sum(log(sigma));

  // Priors (weakly informative)
  mu[1] ~ normal(0, 1);
  mu[2] ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  sigma ~ normal(0, 1); // Half-normal
  rho   ~ normal(0, 0.5);

  // Vectorized calculation of residuals
  // Prediction[t] = mu + Y[t-1] * Phi'
  // rep_matrix broadcasts mu across all T-1 rows.
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals = y[2:T, ] - predictions;

  // Likelihood calculation (Loop required for copula)
  for (t in 1:T-1) {
    row_vector[2] res = residuals[t];
    
    // Optimization: Standardize residuals
    row_vector[2] z = res ./ sigma';

    // 1. Marginal contributions (using optimized standardized form)
    target += std_normal_lpdf(z[1]);
    target += std_normal_lpdf(z[2]);
    // Apply Jacobian adjustment for the scaling transformation
    target += -log_sigma_sum;

    // 2. Copula contribution
    real u1 = Phi(z[1]);
    real u2 = Phi(z[2]);
    target += gaussian_copula_ld(u1, u2, rho);
  }
}

generated quantities {
  // Output the standard Phi matrix for interpretation
  matrix[2, 2] Phi = Phi_T';
}
