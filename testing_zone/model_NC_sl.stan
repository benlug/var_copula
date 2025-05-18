// file: model_NC_sl.stan
// Bivariate VAR(1) model assuming Normal margins and Clayton copula.
// Single-level (no random effects).

functions {
  // Clayton copula log-density for 2 dimensions
  // theta > 0
  real clayton_copula_density_2d(real u, real v, real theta) {
    real log_lik = 0;
    real eps = 1e-9; // Epsilon for clamping

    // Clamp inputs
    real u_clamp = fmax(eps, fmin(1.0 - eps, u));
    real v_clamp = fmax(eps, fmin(1.0 - eps, v));

    if (theta <= eps) return negative_infinity(); // theta must be > 0

    real term1 = (1.0 + theta) * log(theta); // Correction: Stan uses log(theta) * (1+theta) form? No, log(1+theta) is simpler.
    // Simpler form: log(1+theta) - (1+theta)*log(u*v) - (2+1/theta)*log(u^(-theta) + v^(-theta) - 1)
    // Let's use the PDF form C(u,v) = (1+theta) * (u*v)^(-1-theta) * (u^(-theta) + v^(-theta) - 1)^(-2-1/theta)
    // log C(u,v) = log(1+theta) + (-1-theta)*(log(u)+log(v)) + (-2-1/theta)*log(u^(-theta) + v^(-theta) - 1)

    real log_u = log(u_clamp);
    real log_v = log(v_clamp);

    real u_pow_neg_theta = pow(u_clamp, -theta);
    real v_pow_neg_theta = pow(v_clamp, -theta);

    real sum_pow = u_pow_neg_theta + v_pow_neg_theta - 1.0;

    // Check if sum_pow is positive, otherwise density is 0 (log = -inf)
    if (sum_pow <= eps) return negative_infinity();

    log_lik += log1p(theta); // log(1 + theta)
    log_lik += (-1.0 - theta) * (log_u + log_v);
    log_lik += (-2.0 - (1.0 / theta)) * log(sum_pow);

    return log_lik;
  }
}

data {
  int<lower=2> T;       // Number of time points
  vector[2] y[T];       // Data: array of 2-element vectors y_t
}

parameters {
  // Global Parameters Only
  vector[2] mu;                // Intercepts
  real<lower=-1, upper=1> phi11; // VAR param
  real<lower=-1, upper=1> phi12; // VAR param
  real<lower=-1, upper=1> phi21; // VAR param
  real<lower=-1, upper=1> phi22; // VAR param

  // Residual Distribution Parameters (Normal)
  vector<lower=0>[2] sigma;      // Residual SDs (sigma1, sigma2)

  // Copula Parameter (Clayton)
  real<lower=0> theta;           // Clayton copula parameter theta > 0
}

transformed parameters {
  matrix[2, 2] Phi;
  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;
}

model {
  // Priors
  mu ~ normal(0, 1);
  phi11 ~ normal(0, 0.5);
  phi12 ~ normal(0, 0.5);
  phi21 ~ normal(0, 0.5);
  phi22 ~ normal(0, 0.5);
  sigma ~ normal(0, 1); // Half-Normal implied
  theta ~ gamma(2, 1);  // Prior for theta > 0 (adjust shape/rate as needed, gamma(2,1) has mean 2)

  // Likelihood Calculation
  for (t in 2:T) {
    vector[2] y_curr = y[t];
    vector[2] y_prev = y[t-1];
    vector[2] cond_mean = mu + Phi * y_prev;
    vector[2] residuals = y_curr - cond_mean;

    // 1. Add marginal log-likelihoods (Normal)
    target += normal_lpdf(residuals[1] | 0, sigma[1]);
    target += normal_lpdf(residuals[2] | 0, sigma[2]);

    // 2. Calculate CDFs and add Clayton Copula density
    real u1 = normal_cdf(residuals[1] | 0, sigma[1]);
    real u2 = normal_cdf(residuals[2] | 0, sigma[2]);
    target += clayton_copula_density_2d(u1, u2, theta);
  }
}

generated quantities {
  // Optional: Calculate Kendall's Tau from theta
  real<lower=0, upper=1> tau = theta / (theta + 2.0);
}
