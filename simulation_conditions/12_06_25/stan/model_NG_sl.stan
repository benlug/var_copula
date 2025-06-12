functions {
  real gaussian_copula_density(real u, real v, real rho) {
    real eps = 1e-9;
    real uu  = fmax(eps, fmin(1 - eps, u));   // safe copies
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
  vector[2] y[T];
}

parameters {
  vector[2] mu;
  real<lower=-1,upper=1> phi11;
  real<lower=-1,upper=1> phi12;
  real<lower=-1,upper=1> phi21;
  real<lower=-1,upper=1> phi22;
  vector<lower=0>[2] sigma;
  real<lower=-1,upper=1> rho;
}

transformed parameters {
  matrix[2,2] Phi;
  Phi[1,1]=phi11; Phi[1,2]=phi12;
  Phi[2,1]=phi21; Phi[2,2]=phi22;
}

model {
  mu    ~ normal(0,1);
  phi11 ~ normal(0,0.5);
  phi12 ~ normal(0,0.5);
  phi21 ~ normal(0,0.5);
  phi22 ~ normal(0,0.5);
  sigma ~ normal(0,1);
  rho   ~ normal(0,0.5);

  for (t in 2:T) {
    vector[2] res = y[t] - (mu + Phi * y[t-1]);
    target += normal_lpdf(res | 0, sigma);
    real u1 = normal_cdf(res[1] | 0, sigma[1]);
    real u2 = normal_cdf(res[2] | 0, sigma[2]);
    target += gaussian_copula_density(u1, u2, rho);
  }
}
