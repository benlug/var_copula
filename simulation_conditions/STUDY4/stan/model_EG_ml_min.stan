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
  int<lower=1> N;
  int<lower=2> T;
  array[N] matrix[T,2] y;
  vector[2] skew_direction;   // +/-1 per series (shared across units)
}

parameters {
  // Hierarchy on intercepts (minimal pooling)
  row_vector[2] mu_bar;
  vector<lower=0>[2] tau_mu;
  matrix[N,2] z_mu;           // non-centered: mu[i]=mu_bar + tau .* z

  // VAR dynamics (global)
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // EG margins: global slack over feasibility bound
  vector[2] eta;              // sigma_exp = b + exp(eta)

  // Gaussian copula correlation (global)
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2,2] Phi_T;
  matrix[N,2] mu;

  Phi_T[1,1] = phi11;  Phi_T[2,1] = phi12;
  Phi_T[1,2] = phi21;  Phi_T[2,2] = phi22;

  for (i in 1:N) {
    mu[i] = mu_bar + (to_row_vector(tau_mu) .* z_mu[i]);
  }
}

model {
  vector[2] b;           // global feasibility bounds across all (i,t)
  vector[2] sigma_exp_local;
  vector[2] rate_exp;

  // Priors
  mu_bar ~ normal(0, 1);
  to_vector(z_mu) ~ normal(0, 1);
  tau_mu ~ student_t(3, 0, 0.5); // half-t via lower=0 constraint

  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho   ~ normal(0, 0.5);

  // 1) Scan residuals to compute global feasibility bounds b_j
  {
    // Initialize b using first available residual
    row_vector[2] pred0 = mu[1] + y[1][1] * Phi_T;
    row_vector[2] res0  = y[1][2] - pred0;
    for (j in 1:2) b[j] = -skew_direction[j] * res0[j];

    for (i in 1:N) {
      for (t in 2:T) {
        row_vector[2] pred = mu[i] + y[i][t-1] * Phi_T;
        row_vector[2] res  = y[i][t] - pred;
        for (j in 1:2) {
          real v = -skew_direction[j] * res[j];
          if (v > b[j]) b[j] = v;
        }
      }
    }
  }

  // 2) Global EG scales: sigma = b + exp(eta)
  for (j in 1:2) sigma_exp_local[j] = b[j] + exp(eta[j]);
  rate_exp = 1.0 ./ sigma_exp_local;

  // prior on sigma_exp via change of variables: lognormal(0,0.5)
  for (j in 1:2) target += lognormal_lpdf(sigma_exp_local[j] | 0, 0.5) + eta[j];

  // 3) Likelihood over all units and times
  for (i in 1:N) {
    for (t in 2:T) {
      row_vector[2] pred = mu[i] + y[i][t-1] * Phi_T;
      row_vector[2] res  = y[i][t] - pred;

      vector[2] x; vector[2] u;
      for (j in 1:2) {
        x[j] = sigma_exp_local[j] + skew_direction[j] * res[j];  // >0
        target += exponential_lpdf(x[j] | rate_exp[j]);
        u[j] = exponential_cdf(x[j], rate_exp[j]);
        if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
      }
      target += gaussian_copula_ld(u[1], u[2], rho);
    }
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';

  // Recompute and expose sigma_exp for analysis
  // Using consistent naming: sigma_exp[1], sigma_exp[2]
  vector[2] sigma_exp;
  vector[2] b_gq;      // feasibility bound (for diagnostics)
  vector[2] slack_gq;  // exp(eta) slack term (for diagnostics)

  {
    vector[2] b_local;
    
    // Initialize from first residual
    row_vector[2] pred0 = mu[1] + y[1][1] * Phi_T;
    row_vector[2] res0  = y[1][2] - pred0;
    for (j in 1:2) b_local[j] = -skew_direction[j] * res0[j];
    
    // Scan all residuals
    for (i in 1:N) {
      for (t in 2:T) {
        row_vector[2] pred = mu[i] + y[i][t-1] * Phi_T;
        row_vector[2] res  = y[i][t] - pred;
        for (j in 1:2) {
          real v = -skew_direction[j] * res[j];
          if (v > b_local[j]) b_local[j] = v;
        }
      }
    }
    
    for (j in 1:2) {
      b_gq[j]      = b_local[j];
      slack_gq[j]  = exp(eta[j]);
      sigma_exp[j] = b_gq[j] + slack_gq[j];
    }
  }
}
