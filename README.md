# epinow2julia: EpiNow2 powered by Julia/Turing.jl

> **This is an experimental rewrite of
> [EpiNow2](https://github.com/epiforecasts/EpiNow2) that replaces all
> internal Stan models with calls to
> [EpiNow2.jl](https://github.com/epiforecasts/epinow2.jl), a Julia
> implementation using [Turing.jl](https://turinglang.org/) for Bayesian
> inference.** It is not yet released and is under active development.
> See the [feature status](#feature-status) below.

## What is this?

`{epinow2julia}` provides the familiar EpiNow2 interface
(`estimate_infections()`, `epinow()`, `regional_epinow()`) but delegates
all Bayesian inference to Julia/Turing.jl via the
[JuliaConnectoR](https://cran.r-project.org/package=JuliaConnectoR)
package. This eliminates all Stan code from EpiNow2 and uses a pure
Julia implementation for the renewal equation model, observation model,
and MCMC sampling.

## Feature status

| Feature                 | Status  | Notes                                             |
|-------------------------|---------|---------------------------------------------------|
| `estimate_infections()` | Working | Renewal equation via Turing.jl                    |
| `estimate_secondary()`  | Working | Cases to deaths/hospitalisations                  |
| `estimate_truncation()`  | Broken  | Julia-side model issue (AD compatibility)         |
| `simulate_infections()` | Working | Forward simulation from Rt trajectory             |
| `simulate_secondary()`  | Working | Forward simulation of secondary observations      |
| `forecast_infections()` | Working | Forecast from fitted model with new Rt            |
| `epinow()`              | Working | Full pipeline wrapper                             |
| `regional_epinow()`     | Working | Multi-region parallelisation                      |
| `get_samples()`         | Working | R, infections, growth_rate, reported_cases        |
| `get_predictions()`     | Working | summary, sample, quantile formats                 |
| `get_parameters()`      | Working | Returns fitted parameter posteriors                |
| `summary()` / `plot()`  | Working |                                                   |
| Gaussian process on Rt  | Working | `gp_opts()` with Matern/SE kernels               |
| Random walk on Rt       | Working | `rt_opts(rw = 7)` for weekly                     |
| Day-of-week effects     | Working | `obs_opts(week_effect = TRUE)`                    |
| Forecasting             | Working | `forecast_opts(horizon = N)`                      |
| Poisson / NegBin        | Working | `obs_opts(family = ...)`                          |
| Uncertain delays        | Working | Generation time and reporting delays              |
| Population depletion    | Working | `rt_opts(pop = ...)`                              |
| Back-calculation        | Working | `rt_opts(use_rt = FALSE)`                         |
| Stan backend            | Removed | Julia/Turing.jl only                              |

## Installation

This package is not on CRAN. Install from this repository:

```r
# install.packages("pak")
pak::pkg_install("epiforecasts/epinow2julia")
```

### Requirements

- [Julia](https://julialang.org/downloads/) >= 1.10
- The [EpiNow2.jl](https://github.com/epiforecasts/epinow2.jl) Julia
  package (installed in a local project)

### Julia setup

1. Install Julia (e.g. via [juliaup](https://github.com/JuliaLang/juliaup))

2. Clone and set up the Julia package:
   ```bash
   git clone https://github.com/epiforecasts/epinow2.jl.git
   cd epinow2.jl
   julia --project -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
   ```

3. Tell R where to find it:
   ```r
   options(EpiNow2.julia_project = "/path/to/epinow2.jl")
   # or set the EPINOW2_JULIA_PROJECT environment variable
   ```

## Quick start

```r
library(EpiNow2)
options(EpiNow2.julia_project = "/path/to/epinow2.jl")

# Example data
reported_cases <- example_confirmed[1:40]

# Generation time and delays
generation_time <- LogNormal(meanlog = 1.6, sdlog = 0.5, max = 14)
incubation_period <- LogNormal(mean = 5, sd = 3, max = 14)
reporting_delay <- LogNormal(mean = 2, sd = 1, max = 10)

# Estimate Rt
result <- estimate_infections(
  reported_cases,
  generation_time = gt_opts(generation_time),
  delays = delay_opts(incubation_period + reporting_delay),
  rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.1)),
  stan = stan_opts(samples = 1000, warmup = 250, chains = 4)
)

summary(result)
plot(result)
```

## How it works

The R package converts all user inputs (distributions, options, data)
into Julia equivalents via JuliaConnectoR, then calls the corresponding
Julia function. Results are converted back to R data.tables in the same
format as the original EpiNow2 package.

```
R user code
    │
    ▼
EpiNow2 R interface (identical API)
    │
    ▼
R/convert.R (R opts → Julia opts)
    │
    ▼
JuliaConnectoR ←→ Julia process
    │
    ▼
EpiNow2.jl / Turing.jl (MCMC)
    │
    ▼
R/convert.R (Julia results → R data.tables)
    │
    ▼
R user gets familiar output
```

## Key differences from EpiNow2

- **Backend**: Julia/Turing.jl instead of Stan (rstan/cmdstanr)
- **First-run latency**: Julia compiles on first use (~1-2 minutes).
  Subsequent calls in the same session are fast.
- **System dependency**: Requires Julia instead of a C++ toolchain
- **`stan_opts()`**: Still accepted for backward compatibility; maps
  internally to Julia inference options
- **Numerical results**: Will differ slightly from Stan due to different
  MCMC implementations, but should be statistically equivalent

## Citation

If using this package, please cite the original EpiNow2 paper:

> Abbott, S. et al. (2020). "Estimating the time-varying reproduction
> number of SARS-CoV-2 using national and subnational case counts."
> *Wellcome Open Research* 5:112.
> doi:[10.12688/wellcomeopenres.16006.1](https://doi.org/10.12688/wellcomeopenres.16006.1)
