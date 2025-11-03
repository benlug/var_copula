// sem_B_latent_EG.stan
// Study B: skew in state innovations (zeta); measurement error ≡ 0.
// Exponential margins for ζ with Gaussian copula rho.
// No latent "state" parameter; the likelihood is written directly for y_t.
// This is equivalent to a VAR(1) with Exponential innovations (SEM view).

functions {
  real exp_margin_lpdf(real e, real s, int dir) {
    if (dir == 1) {
      real x = e / s + 1; if (x < 0) return negative_infinity();
      return exponential_lpdf(x | 1) - log(s);
    } else {
      real x = 1 - e / s; if (x < 0) return negative_infinity();
      return exponential_lpdf(x | 1) - log(s);
    }
  }
  real exp_margin_cdf(real e, real s, int dir) {
    if (dir == 1) return 1 - exp(-(e / s + 1));
    else          return exp(e / s - 1);
  }
  real gauss_copula_lpdf(real u1, real u2, real rho) {
    vector[2] z; matrix[2,2] S; vector[2] mu0;
    z[1]=inv_Phi(fmin(fmax(u1,1e-9),1-1e-9));
    z[2]=inv_Phi(fmin(fmax(u2,1e-9),1-1e-9));
    S[1,1]=1; S[1,2]=rho; S[2,1]=rho; S[2,2]=1; mu0[1]=0; mu0[2]=0;
    return multi_normal_lpdf(z | mu0, S)
         - normal_lpdf(z[1] | 0,1) - normal_lpdf(z[2] | 0,1);
  }
}

data {
  int<lower=1> T;
  matrix[T,2] y;                          // observed series
  int<lower=-1, upper=1> skew_direction[2];
}

parameters {
  vector[2] mu;
  real phi11; real phi12;
  real phi21; real phi22;
  vector[2] eta;                         // log sigma_exp for ζ
  real<lower=-0.995, upper=0.995> rho;   // copula correlation at state layer
}

transformed parameters {
  vector<lower=0>[2] sigma_exp = exp(eta);
  matrix[2,2] B;
  B[1,1]=phi11; B[1,2]=phi12;
  B[2,1]=phi21; B[2,2]=phi22;
}

model {
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.3);  phi21 ~ normal(0, 0.3);
  eta   ~ normal(0, 0.5);
  rho   ~ normal(0, 0.5);

  // Innovations ζ_t implied by y_t = mu + B y_{t-1} + ζ_t (for t >= 2)
  for (t in 2:T) {
    vector[2] mean_t = mu + B * (y[t-1]') ;
    real z1 = y[t,1] - mean_t[1];
    real z2 = y[t,2] - mean_t[2];

    target += exp_margin_lpdf(z1 | sigma_exp[1], skew_direction[1]);
    target += exp_margin_lpdf(z2 | sigma_exp[2], skew_direction[2]);

    real u1 = exp_margin_cdf(z1 | sigma_exp[1], skew_direction[1]);
    real u2 = exp_margin_cdf(z2 | sigma_exp[2], skew_direction[2]);
    target += gauss_copula_lpdf(u1 | u2, rho);
  }
}
