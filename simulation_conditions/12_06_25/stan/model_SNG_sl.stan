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
  int<lower=2> T;          // series length
  vector[2] y[T];          // observed values
}

parameters {
  //---------------- location / AR parameters ----------------------------
  vector[2] mu;
  real<lower=-1,upper=1> phi11;
  real<lower=-1,upper=1> phi12;
  real<lower=-1,upper=1> phi21;
  real<lower=-1,upper=1> phi22;

  //---------------- skew‑normal innovation parameters -------------------
  vector<lower=0>[2] omega;              // scale (>0)
  vector<lower=-1,upper=1>[2] delta;     // NEW: bounded shape ratio (−1,1)

  //---------------- copula correlation ----------------------------------
  real<lower=-1,upper=1> rho;
}

transformed parameters {
  //---------------- AR matrix -------------------------------------------
  matrix[2,2] Phi;
  Phi[1,1] = phi11;  Phi[1,2] = phi12;
  Phi[2,1] = phi21;  Phi[2,2] = phi22;

  //---------------- original alpha & xi (kept for output compatibility) --
  vector[2] alpha;
  vector[2] xi;

  alpha = delta ./ sqrt(1 - square(delta));            // back‑transform
  xi    = -omega .* (delta * sqrt(2 / pi()));          // keeps E[eps]=0
}

model {
  //---------------- priors (same scale as before, but on delta) ----------
  mu     ~ normal(0, 1);
  phi11  ~ normal(0, 0.5);
  phi12  ~ normal(0, 0.5);
  phi21  ~ normal(0, 0.5);
  phi22  ~ normal(0, 0.5);

  omega  ~ normal(0, 1);                  // weakly informative
  ddelta ~ normal(0, 0.5);                // truncated by bounds −1,1
  rho    ~ normal(0, 0.5);

  //---------------- likelihood -------------------------------------------
  for (t in 2:T) {
    vector[2] res = y[t] - (mu + Phi * y[t-1]);

    // univariate skew‑normal contributions
    target += skew_normal_lpdf(res[1] | xi[1], omega[1], alpha[1]);
    target += skew_normal_lpdf(res[2] | xi[2], omega[2], alpha[2]);

    // copula density contribution
    real u1 = skew_normal_cdf(res[1] | xi[1], omega[1], alpha[1]);
    real u2 = skew_normal_cdf(res[2] | xi[2], omega[2], alpha[2]);
    target += gaussian_copula_density(u1, u2, rho);
  }
}
