// sem_A_indicator_EG.stan
// Study A: skew in measurement layer (epsilon), state innovations Gaussian.
// Exponential margins are standardized; scale estimated via sigma_exp = exp(eta).
// Gaussian copula links the two measurement residuals (per time) via rho.
// Quantiles/CDFs per Study‑4 (standardized Exp). :contentReference[oaicite:6]{index=6}

functions {
  // Standardized exponential margins: E(+) = X - 1, X~Exp(1); E(-) = -E(+)
  real exp_margin_lpdf(real e, real s, int dir) {
    // dir = +1 for E(+), dir = -1 for E(-)
    if (dir == 1) {
      real x = e / s + 1;         // support e >= -s
      if (x < 0) return negative_infinity();
      return exponential_lpdf(x | 1) - log(s);
    } else {
      real x = 1 - e / s;         // support e <= s
      if (x < 0) return negative_infinity();
      return exponential_lpdf(x | 1) - log(s);
    }
  }
  real exp_margin_cdf(real e, real s, int dir) {
    if (dir == 1) {
      return 1 - exp(-(e / s + 1));   // e >= -s
    } else {
      return exp(e / s - 1);          // e <=  s
    }
  }
  // Gaussian copula density in log-space (requires u in (0,1))
  real gauss_copula_lpdf(real u1, real u2, real rho) {
    vector[2] z;
    matrix[2,2] S;
    vector[2] mu0;
    z[1] = inv_Phi(fmin(fmax(u1, 1e-9), 1 - 1e-9));
    z[2] = inv_Phi(fmin(fmax(u2, 1e-9), 1 - 1e-9));
    S[1,1] = 1; S[1,2] = rho;
    S[2,1] = rho; S[2,2] = 1;
    mu0[1] = 0; mu0[2] = 0;
    return multi_normal_lpdf(z | mu0, S)
         - normal_lpdf(z[1] | 0, 1) - normal_lpdf(z[2] | 0, 1);
  }
}

data {
  int<lower=1> T;
  matrix[T,2] y;
  int<lower=-1, upper=1> skew_direction[2]; // +1 right, -1 left
}

parameters {
  vector[2] mu;
  real phi11; real phi12;
  real phi21; real phi22;
  vector[2] eta;                         // log sigma_exp
  real<lower=-0.995, upper=0.995> rho;   // copula correlation at measurement layer
  matrix[T,2] state;                     // latent states
}

transformed parameters {
  vector<lower=0>[2] sigma_exp = exp(eta);
  matrix[2,2] B;
  B[1,1]=phi11; B[1,2]=phi12;
  B[2,1]=phi21; B[2,2]=phi22;
}

model {
  // Priors (weakly informative)
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.3);  phi21 ~ normal(0, 0.3);
  eta   ~ normal(0, 0.5);         // sigma_exp ~ lognormal(0,0.5)
  rho   ~ normal(0, 0.5);

  // Gaussian state innovations
  for (j in 1:2) state[1,j] ~ normal(0, 1);
  for (t in 2:T) {
    vector[2] mean_t = mu + B * (state[t-1]') ;
    for (j in 1:2)
      state[t,j] ~ normal(mean_t[j], 1);
  }

  // Measurement residuals: epsilon_t = y_t - state_t with Exp margins + copula
  for (t in 1:T) {
    real e1 = y[t,1] - state[t,1];
    real e2 = y[t,2] - state[t,2];

    // marginals
    target += exp_margin_lpdf(e1 | sigma_exp[1], skew_direction[1]);
    target += exp_margin_lpdf(e2 | sigma_exp[2], skew_direction[2]);

    // copula on (e1,e2)
    real u1 = exp_margin_cdf(e1 | sigma_exp[1], skew_direction[1]);
    real u2 = exp_margin_cdf(e2 | sigma_exp[2], skew_direction[2]);
    target += gauss_copula_lpdf(u1 | u2, rho);
  }
}
