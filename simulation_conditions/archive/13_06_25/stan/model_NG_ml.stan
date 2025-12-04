functions {
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9;
    real uu  = fmax(eps, fmin(1 - eps, u));
    real vv  = fmax(eps, fmin(1 - eps, v));
    real z1  = inv_Phi(uu);
    real z2  = inv_Phi(vv);
    real rho2 = square(rho);
    return -0.5 * log1m(rho2)
           - 0.5 / (1 - rho2) *
             (square(z1) - 2 * rho * z1 * z2 + square(z2))
           + 0.5 * (square(z1) + square(z2));
  }
}

data {
  int<lower=1> N;
  int<lower=2> T;
  vector[2] y[N, T];
}

parameters {
  //----------------  level‑2  --------------------------
  vector[2] mu0;
  vector<lower=0>[2] sigma_mu;

  vector[4] Phi0;
  vector<lower=0>[4] sigma_phi;

  //----------------  level‑1  --------------------------
  vector[2] mu[N];
  vector[4] Phi_vec[N];

  //----------------  innovation  -----------------------
  vector<lower=0>[2] sigma;
  real<lower=-1,upper=1> rho;
}

model {
  // priors
  mu0        ~ normal(0, 10);
  sigma_mu   ~ normal(0, 1);

  Phi0       ~ normal(0, 5);
  sigma_phi  ~ normal(0, 1);

  sigma      ~ normal(1, 1);
  rho        ~ normal(0, 0.5);

  // random effects
  for (i in 1:N) {
    mu[i]      ~ normal(mu0,  sigma_mu);
    Phi_vec[i] ~ normal(Phi0, sigma_phi);
  }

  // likelihood
  for (i in 1:N) {
    matrix[2,2] Phi_i;
    Phi_i[1,1] = Phi_vec[i,1];
    Phi_i[1,2] = Phi_vec[i,2];
    Phi_i[2,1] = Phi_vec[i,3];
    Phi_i[2,2] = Phi_vec[i,4];

    for (t in 2:T) {
      vector[2] res = y[i, t] - (mu[i] + Phi_i * y[i, t-1]);

      target += normal_lpdf(res | 0, sigma);

      real u1 = normal_cdf(res[1] | 0, sigma[1]);
      real u2 = normal_cdf(res[2] | 0, sigma[2]);
      target += gaussian_copula_density(u1, u2, rho);
    }
  }
}
