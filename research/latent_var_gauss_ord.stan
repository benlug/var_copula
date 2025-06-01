# save this as "latent_var_gauss_ord.stan" and then compile

stan_code <- '
data {
  int<lower=2> T;              // number of time points
  vector[T] Y1;                // observed continuous
  int<lower=1> Y2[T];          // observed ordinal, e.g. 1..K
  int<lower=2> K;              // number of categories for Y2
  vector[K-1] thresh;          // known thresholds for ordinal
}
parameters {
  vector[2] mu;                // latent means
  matrix[2,2] Phi;             // VAR(1) coefficients
  cholesky_factor_cov[2] Lres; // Cholesky factor for residual cov (2x2)
  real<lower=0,upper=1e-1> sigma_meas; // small measurement error for Y1

  // latent states Z[t] = (Z1[t], Z2[t])
  // We separate them out for clarity
  vector[T] z1;
  vector[T] z2;
}
transformed parameters {
  cov_matrix[2] Sigma;
  Sigma = multiply_lower_tri_self_transpose(Lres); // Lres*Lres^T
}
model {
  // Priors (very basic)
  mu ~ normal(0, 2);
  to_vector(Phi) ~ normal(0, 0.5);  // fairly loose
  // enforce stationarity? (omitted for brevity, or we do a prior on spectral radius)
  Lres ~ lkj_corr_cholesky(2);

  // small meas error
  sigma_meas ~ normal(0, 0.1);

  // Likelihood:
  // 1) define the dynamic for Z[t], t=2..T
  //    Z[t] ~ Normal( mu + Phi*(Z[t-1] - mu), Sigma )
  for(t in 2:T) {
    vector[2] mean_t;
    mean_t[1] = mu[1] + Phi[1,1]*(z1[t-1]-mu[1]) + Phi[1,2]*(z2[t-1]-mu[2]);
    mean_t[2] = mu[2] + Phi[2,1]*(z1[t-1]-mu[1]) + Phi[2,2]*(z2[t-1]-mu[2]);

    // multi_normal: we store (z1[t], z2[t]) in a vector
    vector[2] zt;
    zt[1] = z1[t];
    zt[2] = z2[t];
    zt ~ multi_normal_cholesky(mean_t, Lres);
  }

  // 2) tie (z1[t]) to observed Y1[t] with small measurement error
  //    => Y1[t] ~ Normal(z1[t], sigma_meas)
  for(t in 1:T) {
    Y1[t] ~ normal(z1[t], sigma_meas);
  }

  // 3) ordinal measurement for Y2:
  //    if Y2[t] == 1 => z2[t] < thresh[1]
  //    if Y2[t] == 2 => thresh[1] <= z2[t] < thresh[2]
  //    ...
  //    if Y2[t] == K => z2[t] >= thresh[K-1]
  for(t in 1:T) {
    if(Y2[t] == 1) {
      target += normal_lcdf(thresh[1] | z2[t], 1e6) // hack: no scale
                 - normal_lcdf(-1e9 | z2[t], 1e6);   // effectively z2 < thresh[1]
    } else if(Y2[t] == K) {
      target += normal_lcdf(1e9 | z2[t], 1e6)
                 - normal_lcdf(thresh[K-1] | z2[t], 1e6);
    } else {
      // category c => thresh[c-1] <= z2[t] < thresh[c]
      int c = Y2[t];
      target += normal_lcdf(thresh[c] | z2[t], 1e6)
                 - normal_lcdf(thresh[c-1] | z2[t], 1e6);
    }
  }
}
'
