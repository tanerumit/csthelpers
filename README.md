# csthelpers

`csthelpers` is an R package for climate stress testing workflows. It provides tools to
compute scenario weights from climate model ensembles, evaluate and compare weighting
methods, aggregate weighted impacts, and generate climate-surface and radial plots for
analysis and reporting.

The package is aimed at both analysts using climate stress testing outputs and developers
extending the weighting and plotting workflow inside this repository.

## Current Scope

The package currently focuses on five areas:

- GCM weighting: institution-based, genealogy-based, and performance-independence style weights
- Scenario-surface weights: multivariate normal, KDE, and copula-based methods
- Evaluation and diagnostics: method comparison, goodness-of-fit, coverage, and effective sample size summaries
- Impact aggregation: weighted summaries across realizations and scenarios
- Visualization: climate response surfaces, KDE/coplanar diagnostics, basin maps, and radial plots

## Installation

This repository does not currently advertise a CRAN release. Install it locally from the
repo root:

```r
install.packages("devtools")
devtools::install()
```

For interactive development:

```r
devtools::load_all()
```

## Quick Start

### Example: institution-based GCM weights

```r
library(csthelpers)

gcm_data <- data.frame(
  scenario = c("SSP2-4.5", "SSP2-4.5", "SSP2-4.5", "SSP5-8.5"),
  institution = c("CSIRO", "CSIRO", "NOAA-GFDL", "MOHC"),
  model = c("ACCESS-CM2", "ACCESS-ESM1-5", "GFDL-CM4", "UKESM1-0-LL"),
  stringsAsFactors = FALSE
)

weights <- compute_gcm_weights_by_institution(gcm_data)
weights[, c("scenario", "institution", "model", "w_inst")]
```

Within each scenario, the institution-based method gives equal total mass to each
institution, then splits that mass across models and duplicate rows from that institution.

### Example: parametric scenario-surface weights

```r
library(csthelpers)

set.seed(1)
ensemble_data <- data.frame(
  scenario = rep(c("SSP1", "SSP2"), each = 30),
  tavg = c(rnorm(30, 0, 1), rnorm(30, 1, 1)),
  prcp = c(rnorm(30, 0, 1), rnorm(30, 0.5, 1))
)

scenario_grid <- expand.grid(
  tavg = seq(-3, 3, length.out = 40),
  prcp = seq(-3, 3, length.out = 40)
)

w_mvn <- compute_scenario_surface_weights_mvn(
  ensemble_data = ensemble_data,
  scenario_grid = scenario_grid,
  ta_col = "tavg",
  pr_col = "prcp",
  group_col = "scenario",
  normalize = TRUE,
  area_weight = "none",
  verbose = FALSE
)

head(w_mvn)
```

This returns a long-format weight surface over the climate grid, which can then be passed
to evaluation and plotting helpers.

## Main Functions

### Weighting

- `compute_gcm_weights_by_institution()`
- `compute_gcm_weights_by_genealogy()`
- `compute_gcm_weights_bma()`
- `compute_scenario_surface_weights_mvn()`
- `compute_scenario_surface_weights_mvnorm()` for backward compatibility
- `compute_scenario_surface_weights_kde()`
- `compute_scenario_surface_weights_cop()`

### Evaluation and comparison

- `evaluate_scenario_surface_weights()`
- `compare_methods()`
- `copula_goodness_of_fit()`
- `validate_weight_sums()`
- `loo_stability()`
- `stress_region_coverage()`
- `tail_coverage_quantile()`
- `n_active_cells()`
- `effective_n_entropy()`
- `effective_n_hhi()`
- `weight_gini()`
- `weight_max()`
- `weighted_dispersion()`

### Impacts and summaries

- `compute_weighted_impacts()`
- `compute_weighted_ensemble_stats()`
- `summarize_model_weights()`

### Plotting

- `climate_surface_base()`
- `climate_surface_gcm_overlay()`
- `plot_kde()`
- `plot_copula_gof()`
- `radial_plot()`
- `plot_basin_base()`
- `plot_basin_point_feature()`

## Repository Layout

- `R/`: package source
- `tests/testthat/`: unit tests
- `man/`: generated documentation
- `inst/extcode/`: longer runnable examples
- `inst/examples/`: exploratory or usage-oriented scripts
- `data/` and `inst/extdata/`: package/reference data

## Development

Run tests from the repo root:

```r
devtools::test()
```

Regenerate documentation after roxygen changes:

```r
devtools::document()
```

Run linting:

```r
lintr::lint_package()
```

## Notes

- The package metadata in `DESCRIPTION` is still minimal and may be refined as the API stabilizes.
- Some scripts under `inst/examples/` and `inst/extcode/` are exploratory workflows rather than polished end-user tutorials.
- Temporary plots created during development belong in `TEMP/` and should not be kept as durable outputs.
