functions {
  /**
   * Log density of the bivariate Clayton copula.
   * c(u,v) = (1+θ) * (u v)^(-(1+θ)) * (u^{-θ} + v^{-θ} - 1)^(-(2 + 1/θ)),  θ > 0
   */
  real clayton_copula_ld(real u, real v, real theta) {
    real eps = 1e-9;
    real uu  = fmax(eps, fmin(1 - eps, u));
    real vv  = fmax(eps, fmin(1 - eps, v));
    real up  = pow(uu, -theta);
    real vp  = pow(vv, -theta);
    real S   = up + vp - 1;
    return log1p(theta)
         - (1 + theta) * (log(uu) + log(vv))
         - (2 + 1/theta) * log(S);
  }
}

data {
  int<lower=2> T;
  matrix[T, 2] y;
}

parameters {
  row_vector[2] mu;

  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  vector<lower=0>[2] sigma;

  // Clayton parameter (θ > 0). Lower bound avoids numerical issues near 0.
  real<lower=1e-6> theta;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  Phi_T[1, 1] = phi11;  Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21;  Phi_T[2, 2] = phi22;
}

model {
  matrix[T-1, 2] predictions;
  matrix[T-1, 2] residuals;

  // Jacobian for scaling of both margins at each time step
  real log_sigma_sum = sum(log(sigma));

  // Priors (weakly-informative; adjust as needed)
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  sigma ~ normal(0, 1);     // Half-normal around ~1
  theta ~ lognormal(0, 1);  // supports e.g. θ ∈ (0.1, ~7) within ~±2 SDs

  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals   = y[2:T, ] - predictions;

  for (t in 1:(T-1)) {
    row_vector[2] res = residuals[t];
    row_vector[2] z   = res ./ sigma';

    // Marginals: standard normal on standardized residuals
    target += std_normal_lpdf(z[1]);
    target += std_normal_lpdf(z[2]);
    target += -log_sigma_sum;

    // Copula on probability scale
    real u1 = Phi(z[1]);
    real u2 = Phi(z[2]);
    target += clayton_copula_ld(u1, u2, theta);
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';
}
