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
}

parameters {
  // Hierarchy on intercepts (minimal pooling)
  row_vector[2] mu_bar;
  vector<lower=0>[2] tau_mu;
  matrix[N,2] z_mu;

  // VAR dynamics (global)
  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // Gaussian margins (global)
  vector<lower=0>[2] sigma;

  // Gaussian copula correlation (global)
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2,2] Phi_T;
  matrix[N,2] mu;
  Phi_T[1,1] = phi11;  Phi_T[2,1] = phi12;
  Phi_T[1,2] = phi21;  Phi_T[2,2] = phi22;
  for (i in 1:N) mu[i] = mu_bar + (to_row_vector(tau_mu) .* z_mu[i]);
}

model {
  real log_sigma_sum = log(sigma[1]) + log(sigma[2]);

  // Priors
  mu_bar ~ normal(0, 1);
  to_vector(z_mu) ~ normal(0, 1);
  tau_mu ~ student_t(3, 0, 0.5);

  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho   ~ normal(0, 0.5);

  sigma ~ normal(0, 1); // Half-normal via lower=0

  // Likelihood
  for (i in 1:N) {
    for (t in 2:T) {
      row_vector[2] pred = mu[i] + y[i][t-1] * Phi_T;
      row_vector[2] res  = y[i][t] - pred;
      row_vector[2] z    = res ./ sigma';

      target += std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum;

      real u1 = Phi(z[1]);
      real u2 = Phi(z[2]);
      target += gaussian_copula_ld(u1, u2, rho);
    }
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';
}
