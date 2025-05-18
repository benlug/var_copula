# var_copula

## Simulation purpose

- The project focuses on intensive longitudinal data modeled via a multilevel VAR(1) process. The documentation states that the simulation "aims to focus on *skew normal margins* in the context of a multilevel VAR(1)" with the goal of quantifying the impact of ignoring skewness in the innovations.

## Data-generating process (DGP)

- The innovations are drawn from a Gaussian copula with skew normal margins. Each time point uses:
  1. Draw \(\mathbf{z}_{it}\) from a multivariate normal with correlation matrix \(\mathbf{R}\).
  2. Transform each element via the inverse skew normal CDF with shape \(\alpha_j\) to obtain \(\varepsilon_{it,j}\).
  3. Consider left-skew (\(\alpha_j<0\)) and right-skew (\(\alpha_j>0\)) configurations.
  These innovations drive the VAR(1) process, and residual variance/shape parameters are chosen so variability is comparable across conditions.

- The simulation factors include skewness direction, series length \(T\), number of subjects \(N\), autoregressive and cross-lag coefficients, inclusion of random means/lags, and residual variance adjustments to achieve comparable autocorrelation across conditions.

- In code, residuals are produced by sampling from either a Gaussian or Clayton copula and transforming margins via normal or skew-normal quantile functions. The function `generate_residuals` implements this mechanism. These residuals are inserted into a bivariate VAR(1) recursion to generate the time series for each replication.

- The `run_pipeline.R` script defines factor levels: both Gaussian and Clayton copulas, skewness parameters \(\alpha_1,\alpha_2\) at \(-5\) or \(5\), Kendall's \(\tau\) values 0.2 and 0.5, series lengths 30 or 100, AR coefficients \(\phi_{11}=\phi_{22}=0.6\), cross-lags \(\phi_{12}=\phi_{21}=0.2\), and 24 replications per condition.

## Models fitted

- Two classes of models are compared. The documentation lists the fitted models as:
  1. A "Naive Normal-Multilevel-Var(1)" assuming normal innovations.
  2. A "Copula-Based Model" correctly specifying skew normal margins and matching the DGP's copula.

- Corresponding Stan programs implement four single-level variants: Normal margins with Gaussian copula (NG), Normal with Clayton (NC), Skew-Normal with Gaussian (SG), and Skew-Normal with Clayton (SC). `fit_models.R` fits the correct specification for each dataset and the standard normal–Gaussian baseline.

## Goal of the study

- The summary section clarifies the aim: the simulation "systematically varies skewness, autocorrelation, sample size, and random effects," comparing a naive normal-based DSEM with the correctly specified copula-based model to evaluate how misspecification affects bias, coverage, Type I error, and diagnostics. Metrics include bias, RMSE, coverage rates, convergence checks, and power for cross-lag effects.

In essence, the study generates bivariate VAR(1) data where innovations follow skew normal margins coupled via either a Gaussian or Clayton copula. It then fits both standard normal-based and copula-based models to quantify the consequences of ignoring skewness and to assess estimation performance across various experimental factors.
