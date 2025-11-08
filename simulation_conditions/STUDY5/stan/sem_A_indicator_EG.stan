// sem_A_indicator_EG.stan  (NON-CENTERED WITH FIX 3)
// Study A: indicator-skew (skew in measurement errors), Gaussian state dynamics.
// Exponential margins for measurement residuals with Gaussian copula at the
// active layer. Non-centered state innovations to improve mixing on long T.
// FIX 3: Added bounds checking to detect and reject support violations
//
// References: Study‑4 brief for model placement; Study‑2 coding style. 

functions {
  // One‑sided standardized exponential margin (mean 0, sd 1) with scale s>0.
  // dir = +1 (right-skew): support e >= -s ; dir = -1 (left-skew): support e <= s.
  real exp_margin_lpdf(real e, real s, int dir) {
    if (dir == 1) {
      real x = e / s + 1;             // x >= 0 required
      if (x < 0) return negative_infinity();
      return exponential_lpdf(x | 1) - log(s);  // Jacobian
    } else {
      real x = 1 - e / s;             // x >= 0 required
      if (x < 0) return negative_infinity();
      return exponential_lpdf(x | 1) - log(s);  // Jacobian
    }
  }
  real exp_margin_cdf(real e, real s, int dir) {
    // safe, monotone CDFs (no log terms)
    if (dir == 1) return 1 - exp(-(e / s + 1));
    else          return exp(e / s - 1);
  }
  // Gaussian copula log-density for scalar uniforms u1,u2 in (0,1).
  real gauss_copula_lpdf(real u1, real u2, real rho) {
    vector[2] z; matrix[2,2] S; vector[2] mu0;
    z[1] = inv_Phi(fmin(fmax(u1, 1e-9), 1 - 1e-9));
    z[2] = inv_Phi(fmin(fmax(u2, 1e-9), 1 - 1e-9));
    S[1,1] = 1; S[1,2] = rho; S[2,1] = rho; S[2,2] = 1;
    mu0[1] = 0; mu0[2] = 0;
    return multi_normal_lpdf(z | mu0, S)
         - normal_lpdf(z[1] | 0, 1) - normal_lpdf(z[2] | 0, 1);
  }
}

data {
  int<lower=1> T;               // time points
  matrix[T,2] y;                // observed indicators
  int<lower=-1, upper=1> skew_direction[2];  // +1 right, -1 left per margin
}

parameters {
  // VAR(1) on the latent state (Gaussian innovations)
  vector[2] mu;
  real phi11; real phi12;
  real phi21; real phi22;

  // measurement Exponential scales (log scale)
  vector[2] eta;   // sigma_exp = exp(eta)

  // same-time copula correlation at measurement layer
  real<lower=-0.995, upper=0.995> rho;

  // NON-CENTERED innovations for the latent state evolution
  // z_raw[t,] ~ N(0, I), and we build state deterministically from these.
  matrix[T,2] z_raw;
}

transformed parameters {
  vector<lower=0>[2] sigma_exp = exp(eta);
  matrix[2,2] B;
  matrix[T,2] state;

  B[1,1] = phi11; B[1,2] = phi12;
  B[2,1] = phi21; B[2,2] = phi22;

  {
    vector[2] s;
    // state recursion using non-centered shocks
    for (t in 1:T) {
      vector[2] z_t = to_vector(z_raw[t]);   // standard normal shocks
      if (t == 1) {
        s = mu + z_t;                        // prior: state[1] ~ N(mu, I)
      } else {
        s = mu + B * s + z_t;                // state[t] = mu + B*state[t-1] + z_t
      }
      state[t,1] = s[1];
      state[t,2] = s[2];
    }
  }
}

model {
  // priors (weakly informative)
  mu    ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);  phi22 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.3);  phi21 ~ normal(0, 0.3);

  // keep measurement scales from collapsing too small at init;
  // still wide enough for data to adjust.
  eta   ~ normal(log(1.0), 0.7);

  rho   ~ normal(0, 0.5);

  // non-centered state innovations
  // FIX 3: Explicit prior on z_raw to prevent extreme values
  to_vector(z_raw) ~ normal(0, 1);

  // measurement likelihood with Gaussian copula at the active layer
  for (t in 1:T) {
    real e1 = y[t,1] - state[t,1];
    real e2 = y[t,2] - state[t,2];
    
    // FIX 3: Check support violations before computing likelihood
    // This provides better diagnostics than silent failures
    if (skew_direction[1] == 1) {
      if (e1 < -sigma_exp[1]) {
        reject("Support violation at t=", t, ": e1=", e1, 
               " < -sigma_exp[1]=", -sigma_exp[1], 
               " (right-skew requires e1 >= -sigma_exp[1])");
      }
    } else {
      if (e1 > sigma_exp[1]) {
        reject("Support violation at t=", t, ": e1=", e1, 
               " > sigma_exp[1]=", sigma_exp[1], 
               " (left-skew requires e1 <= sigma_exp[1])");
      }
    }
    
    if (skew_direction[2] == 1) {
      if (e2 < -sigma_exp[2]) {
        reject("Support violation at t=", t, ": e2=", e2, 
               " < -sigma_exp[2]=", -sigma_exp[2], 
               " (right-skew requires e2 >= -sigma_exp[2])");
      }
    } else {
      if (e2 > sigma_exp[2]) {
        reject("Support violation at t=", t, ": e2=", e2, 
               " > sigma_exp[2]=", sigma_exp[2], 
               " (left-skew requires e2 <= sigma_exp[2])");
      }
    }
    
    // marginal terms (one-sided support enforced by exp_margin_lpdf)
    target += exp_margin_lpdf(e1 | sigma_exp[1], skew_direction[1]);
    target += exp_margin_lpdf(e2 | sigma_exp[2], skew_direction[2]);

    // copula term
    {
      real u1 = exp_margin_cdf(e1 | sigma_exp[1], skew_direction[1]);
      real u2 = exp_margin_cdf(e2 | sigma_exp[2], skew_direction[2]);
      target += gauss_copula_lpdf(u1 | u2, rho);
    }
  }
}
