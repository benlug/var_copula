// model_nc_sl.stan
// single-level bivariate var(1) with normal margins
// dependence between residuals is modeled via a clayton copula

functions {
  real clayton_copula_density_2d(real u, real v, real theta) {
    real log_lik = 0;
    real eps = 1e-9; // epsilon for clamping

    // clamp inputs
    real u_clamp = fmax(eps, fmin(1.0 - eps, u));
    real v_clamp = fmax(eps, fmin(1.0 - eps, v));

    if (theta <= eps) return negative_infinity(); // theta must be > 0

    real log_u = log(u_clamp);
    real log_v = log(v_clamp);

    real u_pow_neg_theta = pow(u_clamp, -theta);
    real v_pow_neg_theta = pow(v_clamp, -theta);

    real sum_pow = u_pow_neg_theta + v_pow_neg_theta - 1.0;

    // check if sum_pow is positive, otherwise density is 0 (log = -inf)
    if (sum_pow <= eps) return negative_infinity();

    log_lik += log1p(theta); // log(1 + theta)
    log_lik += (-1.0 - theta) * (log_u + log_v);
    log_lik += (-2.0 - (1.0 / theta)) * log(sum_pow);

    return log_lik;
  }
}

data {
  int<lower=2> T;       // number of time points
  vector[2] y[T];       // observed series y[t]
}

parameters {
  vector[2] mu;                // intercepts for y1 and y2
  real<lower=-1, upper=1> phi11; // effect of y1[t-1] on y1[t]
  real<lower=-1, upper=1> phi12; // effect of y2[t-1] on y1[t]
  real<lower=-1, upper=1> phi21; // effect of y1[t-1] on y2[t]
  real<lower=-1, upper=1> phi22; // effect of y2[t-1] on y2[t]

  // residual standard deviations for the normal margins
  vector<lower=0>[2] sigma;

  // clayton copula dependence parameter
  real<lower=0> theta;
}

transformed parameters {
  matrix[2, 2] Phi;
  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;
}

model {
  // priors
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  sigma ~ normal(0, 1); // half-normal implied
  theta ~ gamma(2, 1);  // prior for theta > 0 (gamma mean 2)

  // likelihood
  for (t in 2:T) {
    vector[2] y_curr = y[t];
    vector[2] y_prev = y[t-1];
    vector[2] cond_mean = mu + Phi * y_prev;
    vector[2] residuals = y_curr - cond_mean;

    // add marginal log-likelihoods (normal)
    target += normal_lpdf(residuals[1] | 0, sigma[1]);
    target += normal_lpdf(residuals[2] | 0, sigma[2]);

    // calculate cdfs and add clayton copula density
    real u1 = normal_cdf(residuals[1] | 0, sigma[1]);
    real u2 = normal_cdf(residuals[2] | 0, sigma[2]);
    target += clayton_copula_density_2d(u1, u2, theta);
  }
}

generated quantities {
  real<lower=0, upper=1> tau = theta / (theta + 2.0);
}
