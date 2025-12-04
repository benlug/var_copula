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
  int<lower=1> N;             // persons
  int<lower=2> T;             // occasions
  vector[2] y[N, T];          // series (dim: N × T × 2)
}

parameters {
  //----------------  level‑2  (hyper) -----------------
  vector[2] mu0;
  vector<lower=0>[2] sigma_mu;

  vector[4] Phi0;             // phi11, phi12, phi21, phi22
  vector<lower=0>[4] sigma_phi;

  //----------------  level‑1  (subject‑specific) ------
  vector[2] mu[N];
  vector[4] Phi_vec[N];       // flattened Φ_i

  //----------------  innovation distribution ----------
  vector<lower=0>[2] omega;
  vector[2] alpha;
  real<lower=-1,upper=1> rho;
}

transformed parameters {
  vector[2] xi;
  xi = -omega .* (alpha ./ sqrt(1 + square(alpha))) * sqrt(2 / pi());
}

model {
  // ------------------ priors -------------------------
  mu0        ~ normal(0, 10);
  sigma_mu   ~ normal(0, 1);

  Phi0       ~ normal(0, 5);
  sigma_phi  ~ normal(0, 1);

  omega      ~ normal(1, 1);
  alpha      ~ cauchy(0, 5);
  rho        ~ normal(0, 0.5);

  // ----------------- random effects ------------------
  for (i in 1:N) {
    mu[i]      ~ normal(mu0,  sigma_mu);
    Phi_vec[i] ~ normal(Phi0, sigma_phi);
  }

  // ----------------- likelihood ----------------------
  for (i in 1:N) {
    matrix[2,2] Phi_i;
    Phi_i[1,1] = Phi_vec[i,1];
    Phi_i[1,2] = Phi_vec[i,2];
    Phi_i[2,1] = Phi_vec[i,3];
    Phi_i[2,2] = Phi_vec[i,4];

    for (t in 2:T) {
      vector[2] res = y[i, t] - (mu[i] + Phi_i * y[i, t-1]);

      // univariate skew‑normal contributions
      target += skew_normal_lpdf(res[1] | xi[1], omega[1], alpha[1]);
      target += skew_normal_lpdf(res[2] | xi[2], omega[2], alpha[2]);

      // copula density
      real u1 = skew_normal_cdf(res[1] | xi[1], omega[1], alpha[1]);
      real u2 = skew_normal_cdf(res[2] | xi[2], omega[2], alpha[2]);
      target += gaussian_copula_density(u1, u2, rho);
    }
  }
}
