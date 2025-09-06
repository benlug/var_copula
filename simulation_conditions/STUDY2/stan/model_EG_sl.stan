// model_EG_sl.stan (FIXED: With more robust prior)
// VAR(1) with Shifted Exponential Marginals (estimating scale) and Gaussian Copula.

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
  vector[2] skew_direction; 
}

parameters {
  row_vector[2] mu;
  
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;
  
  vector<lower=0>[2] sigma_exp; 
  
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  vector[2] rate_exp = 1.0 ./ sigma_exp;

  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;
}

model {
  matrix[T-1, 2] residuals;
  matrix[T-1, 2] predictions;
  
  // Priors (weakly informative)
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  
  // ** KEY CHANGE: Use a strictly positive prior for the scale parameter **
  // This prior centers sigma_exp around 1 but prevents it from approaching 0.
  sigma_exp ~ lognormal(0, 0.5); 

  rho ~ normal(0, 0.5);

  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals = y[2:T, ] - predictions;

  for (t in 1:T-1) {
    row_vector[2] res = residuals[t];
    vector[2] x_shifted;
    vector[2] u_vec;
    
    for (i in 1:2) {
      x_shifted[i] = sigma_exp[i] + skew_direction[i] * res[i];
      
      // The intelligent initialization makes this check much less likely to fail,
      // but it remains as a safeguard.
      if (x_shifted[i] <= 0) {
        target += negative_infinity();
        break; 
      }
      
      target += exponential_lpdf(x_shifted[i] | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted[i], rate_exp[i]);
      
      if (skew_direction[i] < 0) {
        u_vec[i] = 1.0 - u_vec[i];
      }
    }
    
    if (min(x_shifted) > 0) {
      target += gaussian_copula_ld(u_vec[1], u_vec[2], rho);
    }
  }
}

generated quantities {
  matrix[2, 2] Phi = Phi_T';
}
