functions {
  real gaussian_copula_ld(real u, real v, real rho) {
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
  matrix[T, 2] y;
  vector[2] skew_direction;      // +1 for right skew, -1 for left skew
  real<lower=0> shape_gam;       // fixed Gamma shape parameter (k)
}

parameters {
  row_vector[2] mu;

  real<lower=-1, upper=1> phi11;
  real<lower=-1, upper=1> phi12;
  real<lower=-1, upper=1> phi21;
  real<lower=-1, upper=1> phi22;

  // Unconstrained parameter for reparameterized residual SD
  // sigma_gam = b + exp(eta), where b is the data-dependent feasibility bound
  vector[2] eta;

  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  Phi_T[1,1] = phi11;  Phi_T[2,1] = phi12;
  Phi_T[1,2] = phi21;  Phi_T[2,2] = phi22;
}

model {
  matrix[T-1, 2] predictions;
  matrix[T-1, 2] residuals;
  vector[2] b;           // feasibility bound for sigma_gam
  vector[2] sigma_gam;   // implied residual SD = b + exp(eta)
  vector[2] rate_gam;    // Gamma rate per margin

  real sqrt_shape = sqrt(shape_gam);

  // Priors
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  rho   ~ normal(0, 0.5);

  // IMPORTANT: As in model_EG_sl.stan, we do NOT place an additional prior on eta.
  // We induce the intended prior on sigma_gam via change-of-variables below.

  // Compute residuals
  predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
  residuals   = y[2:T, ] - predictions;

  // Feasibility bound:
  // For each margin i, define x_shifted = mean_x + s_i * res.
  // We model x_shifted ~ Gamma(shape_gam, rate_gam), where
  //   mean_x = sqrt(shape_gam) * sigma_gam,
  //   var_x  = sigma_gam^2.
  // To ensure x_shifted > 0 for all t we need:
  //   sqrt(shape_gam) * sigma_gam > max_t(-s_i * res_{t,i}).
  // Hence the bound on sigma_gam is:
  //   b_i = max_t(-s_i * res_{t,i}) / sqrt(shape_gam).
  for (i in 1:2) {
    real m = -skew_direction[i] * residuals[1, i];
    for (t in 2:(T-1)) {
      m = fmax(m, -skew_direction[i] * residuals[t, i]);
    }
    b[i] = m / sqrt_shape;
  }

  // Reparameterized sigma_gam: sigma_gam = b + exp(eta) (always feasible)
  for (i in 1:2) {
    sigma_gam[i] = b[i] + exp(eta[i]);
  }

  // Gamma rate so that Var(x_shifted)=sigma_gam^2 and E[x_shifted]=sqrt(shape)*sigma_gam
  rate_gam = sqrt_shape ./ sigma_gam;

  // Induced prior on sigma_gam via change-of-variables
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_gam[i] | 0, 0.5) + eta[i];
  }

  // Likelihood
  for (t in 1:(T-1)) {
    row_vector[2] res = residuals[t];
    vector[2] x_shifted;
    vector[2] u_vec;

    for (i in 1:2) {
      real mean_x = sqrt_shape * sigma_gam[i];
      x_shifted[i] = mean_x + skew_direction[i] * res[i];

      target += gamma_lpdf(x_shifted[i] | shape_gam, rate_gam[i]);
      u_vec[i] = gamma_cdf(x_shifted[i], shape_gam, rate_gam[i]);

      // For left-skewed (mirrored) margins, flip the CDF
      if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
    }

    target += gaussian_copula_ld(u_vec[1], u_vec[2], rho);
  }
}

generated quantities {
  matrix[2,2] Phi = Phi_T';

  vector[2] sigma_gam;
  vector[2] b_gq;
  vector[2] slack;
  vector[2] rate_gam;

  {
    real sqrt_shape = sqrt(shape_gam);
    matrix[T-1, 2] predictions = rep_matrix(mu, T-1) + y[1:T-1, ] * Phi_T;
    matrix[T-1, 2] residuals   = y[2:T, ] - predictions;

    for (i in 1:2) {
      real m = -skew_direction[i] * residuals[1, i];
      for (t in 2:(T-1)) m = fmax(m, -skew_direction[i] * residuals[t, i]);
      b_gq[i]      = m / sqrt_shape;
      slack[i]     = exp(eta[i]);
      sigma_gam[i] = b_gq[i] + slack[i];
      rate_gam[i]  = sqrt_shape / sigma_gam[i];
    }
  }
}
