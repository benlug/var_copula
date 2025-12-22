# Study 6: Time-Varying Parameter Copula-VAR

## Detecting Dynamic Dependence in Intensive Longitudinal Data

### Overview

This simulation study investigates **time-varying copula parameters** in bivariate VAR(1) models for psychological intensive longitudinal data. The key innovation is a state-space formulation where the copula correlation ρ can evolve over time.

### Research Questions

1. **Detection Power**: Under what conditions (T, Δρ, pattern type) can we reliably detect time-varying ρ?
2. **Bias from Ignoring TVP**: How biased are constant-ρ estimates when ρ actually varies?
3. **False Positive Rate**: When ρ is truly constant, how often do TVP models suggest variation?
4. **Computational Feasibility**: Can these models run on typical ESM sample sizes (T = 50-200)?

### Statistical Framework

#### State-Space Model for Time-Varying ρ

```
State equation:    z_t = z_{t-1} + η_t,    η_t ~ N(0, σ_z²)
Observation link:  ρ_t = tanh(z_t)
```

The key parameter `σ_z` controls how much ρ can change:
- σ_z ≈ 0: Constant coupling (no TVP)
- σ_z > 0: Time-varying coupling

#### Simulated Patterns

| Pattern | Description | Psychological Example |
|---------|-------------|----------------------|
| Constant | ρ_t = 0.5 | Stable baseline coupling |
| Linear | ρ: 0.7 → 0.2 | Therapy reducing symptom coupling |
| Step | ρ: 0.3 → 0.7 at midpoint | Episode onset |
| Random Walk | σ_z = 0.05 | General nonstationarity |

### Directory Structure

```
STUDY6_TVP/
├── stan/
│   ├── model_TVP_NG_sl.stan      # Time-varying ρ, Normal margins
│   ├── model_TVP_EG_sl.stan      # Time-varying ρ, Exponential margins
│   ├── model_Const_NG_sl.stan    # Constant ρ, Normal margins (baseline)
│   └── model_Const_EG_sl.stan    # Constant ρ, Exponential margins
├── data/                          # Simulated datasets
├── fits/                          # Stan fit objects
├── results/                       # Analysis outputs
├── simulate_data_tvp.R           # Data generation
├── fit_models_tvp.R              # Model fitting
├── run_pipeline_tvp.R            # Main pipeline
├── analysis_tvp.R                # Result analysis
├── study_6.qmd                   # Quarto report
└── README.md
```

### Running the Pipeline

```r
# Full pipeline
source("run_pipeline_tvp.R")

# Or step by step:
# 1. Generate data
source("simulate_data_tvp.R")
# 2. Fit models
source("fit_models_tvp.R")
# 3. Analyze results
source("analysis_tvp.R")

# Generate report
quarto::quarto_render("study_6.qmd")
```

### Design Factors

| Factor | Levels |
|--------|--------|
| T (time points) | 100, 200 |
| TVP Pattern | Constant, Linear, Step, Random Walk |
| Marginal Distribution | Normal, Exponential |
| VAR Structure | Symmetric (φ = 0.4, cross-lag = 0.1) |
| Replications | 200 per cell |

### Key Output Files

- `summary_replications.csv`: Parameter estimates per replication
- `summary_conditions.csv`: Aggregated metrics per condition
- `summary_sigma_z.csv`: TVP detection analysis
- `summary_rho_recovery.csv`: ρ trajectory recovery metrics

### Evaluation Metrics

#### TVP Detection
- **σ_z posterior**: Does 95% CI exclude 0?
- **Bayes factor**: TVP vs. Constant model (via LOO-CV)

#### Parameter Recovery
- **Bias**: (estimate - truth)
- **Coverage**: 95% CI contains truth?
- **RMSE**: Root mean squared error

#### ρ Trajectory Recovery (TVP models only)
- **ρ_mean recovery**: Matches true mean(ρ_t)?
- **ρ_range recovery**: Captures true variation?

### Dependencies

```r
install.packages(c(
  "rstan",
  "dplyr", "tidyr", "readr", "stringr",
  "ggplot2", "patchwork",
  "doParallel", "foreach",
  "digest"
))
```

### Computational Notes

- **adapt_delta**: 0.9
- **max_treedepth**: 12
- **Runtime**: ~5-10x longer than constant-ρ models
- **Memory**: T-dimensional state vector increases memory requirements
- **Storage**: Fits saved as compact summary objects (~10-50 KB each vs ~5-50 MB for full stanfit)

### Expected Findings

1. **Detection threshold**: TVP detectable when Δρ ≥ 0.3 and T ≥ 100
2. **Pattern effects**: Step changes easier to detect than gradual drift
3. **Φ robustness**: VAR dynamics largely unaffected by ρ misspecification
4. **False positive control**: Shrinkage prior keeps FP rate < 5%

### Citation

Part of the Copula-VAR simulation study program for psychological methodology.

### Contact

[Your contact information]
