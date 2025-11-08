// sem_B_latent_EG.stan  (UPDATED WITH FIXES)
// Study B: latent-skew (skew in state innovations), measurement error ≡ 0.
// FIX: Better support violation handling and more efficient vector operations

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
  matrix[T,2] y;
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

  // y_t = mu + B y_{t-1} + ζ_t   with ζ having exponential margins
  for (t in 2:T) {
    // FIX: More efficient - use to_vector instead of transpose
    vector[2] mean_t = mu + B * to_vector(y[t-1]);
    real z1 = y[t,1] - mean_t[1];
    real z2 = y[t,2] - mean_t[2];

    // FIX: Better support checking with proper error reporting
    // Check support violations and reject with informative message
    if (skew_direction[1] == 1) {
      if (z1 < -sigma_exp[1]) {
        reject("Support violation at t=", t, ": z1=", z1, 
               " < -sigma_exp[1]=", -sigma_exp[1], 
               " (right-skew requires z1 >= -sigma_exp[1])");
      }
    } else {
      if (z1 > sigma_exp[1]) {
        reject("Support violation at t=", t, ": z1=", z1, 
               " > sigma_exp[1]=", sigma_exp[1], 
               " (left-skew requires z1 <= sigma_exp[1])");
      }
    }
    
    if (skew_direction[2] == 1) {
      if (z2 < -sigma_exp[2]) {
        reject("Support violation at t=", t, ": z2=", z2, 
               " < -sigma_exp[2]=", -sigma_exp[2], 
               " (right-skew requires z2 >= -sigma_exp[2])");
      }
    } else {
      if (z2 > sigma_exp[2]) {
        reject("Support violation at t=", t, ": z2=", z2, 
               " > sigma_exp[2]=", sigma_exp[2], 
               " (left-skew requires z2 <= sigma_exp[2])");
      }
    }

    target += exp_margin_lpdf(z1 | sigma_exp[1], skew_direction[1]);
    target += exp_margin_lpdf(z2 | sigma_exp[2], skew_direction[2]);

    {
      // FIX: Compute CDFs only if on-support (which we've now verified above)
      // This is now guaranteed to be valid after the support checks
      real u1 = exp_margin_cdf(z1 | sigma_exp[1], skew_direction[1]);
      real u2 = exp_margin_cdf(z2 | sigma_exp[2], skew_direction[2]);
      target += gauss_copula_lpdf(u1 | u2, rho);
    }
  }
}
