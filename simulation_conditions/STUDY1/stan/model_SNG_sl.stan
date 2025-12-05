functions {
  // Gaussian Copula Log Density
  // Renamed from _lpdf to _ld
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
  matrix[T, 2] y; // Input format updated to matrix
}

transformed data {
  // Define constant for centering calculation: sqrt(2/pi).
  // Using Stan's built-in pi() function is preferred over hardcoded floats.
  real SQRT_2_OVER_PI = sqrt(2.0 / pi()); 
}

parameters {
  // Defined as row_vector for vectorized operations
  row_vector[2] mu;
  
  // AR parameters
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // Skew-normal parameters (Centered Parameterization approach)
  vector<lower=0>[2] omega; // scale (>0)
  // Bounded shape ratio (-1,1) provides stable sampling
  vector<lower=-1, upper=1>[2] delta;     

  // Copula correlation
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2, 2] Phi_T; // Transposed Phi matrix
  vector[2] alpha;    // DP Shape
  vector[2] xi;       // DP Location

  // Construct the transpose
  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;

  // Transform from CP (delta, omega, mu=0) to DP (xi, omega, alpha) 
  // required by Stan's skew_normal_* functions.
  
  // 1. Calculate alpha (shape)
  // alpha = delta / sqrt(1 - delta^2)
  alpha = delta ./ sqrt(1 - square(delta));
  
  // 2. Calculate xi (location) such that E[innovation] = 0
  // xi = -omega * delta * sqrt(2/pi)
  xi = -omega .* (delta * SQRT_2_OVER_PI);
}

model {
  matrix[T-1, 2] residuals;
  matrix[T-1, 2] predictions;

  // Priors
  mu[1] ~ normal(0, 1);
  mu[2] ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);

  omega ~ normal(0, 1);   // Half-normal (omega > 0)
  
  // Prior on delta: regularizes toward symmetry while allowing substantial skewness
  // This prior places more mass near 0 (symmetric) but permits values across (-1, 1)
  delta ~ normal(0, 0.5);
  
  rho ~ normal(0, 0.5);

  // Vectorized calculation of residuals
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals = y[2:T, ] - predictions;

  // Likelihood calculation
  for (t in 1:T-1) {
    row_vector[2] res = residuals[t];

    // 1. Marginal contributions (Skew-Normal PDF)
    target += skew_normal_lpdf(res[1] | xi[1], omega[1], alpha[1]);
    target += skew_normal_lpdf(res[2] | xi[2], omega[2], alpha[2]);

    // 2. Copula contribution (Calculate CDFs)
    real u1 = skew_normal_cdf(res[1], xi[1], omega[1], alpha[1]);
    real u2 = skew_normal_cdf(res[2], xi[2], omega[2], alpha[2]);
    
    // 3. Add copula log density
    target += gaussian_copula_ld(u1, u2, rho);
  }
}

generated quantities {
  // Output the standard Phi matrix for interpretation
  matrix[2, 2] Phi = Phi_T';
}
