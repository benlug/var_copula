functions {
  real gaussian_copula_density(real u, real v, real rho) {
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
  vector[2] y[T];
}

parameters {
  vector[2] mu;
  real<lower=-1,upper=1> phi11;
  real<lower=-1,upper=1> phi12;
  real<lower=-1,upper=1> phi21;
  real<lower=-1,upper=1> phi22;
  vector<lower=0>[2] omega;
  vector[2] alpha;
  real<lower=-1,upper=1> rho;
}

transformed parameters {
  matrix[2,2] Phi;
  vector[2] xi;
  Phi[1,1]=phi11; Phi[1,2]=phi12;
  Phi[2,1]=phi21; Phi[2,2]=phi22;
  xi = -omega .* (alpha ./ sqrt(1 + square(alpha))) * sqrt(2 / pi());
}

model {
  mu    ~ normal(0,1);
  phi11 ~ normal(0,0.5);
  phi12 ~ normal(0,0.5);
  phi21 ~ normal(0,0.5);
  phi22 ~ normal(0,0.5);
  omega ~ normal(0,1);
  alpha ~ cauchy(0,5);
  rho   ~ normal(0,0.5);

  for (t in 2:T) {
    vector[2] res = y[t] - (mu + Phi * y[t-1]);
    target += skew_normal_lpdf(res[1] | xi[1], omega[1], alpha[1]);
    target += skew_normal_lpdf(res[2] | xi[2], omega[2], alpha[2]);
    real u1 = skew_normal_cdf(res[1] | xi[1], omega[1], alpha[1]);
    real u2 = skew_normal_cdf(res[2] | xi[2], omega[2], alpha[2]);
    target += gaussian_copula_density(u1, u2, rho);
  }
}
