


# ------------------------------------------------------------------------------
# 0) Libraries
# ------------------------------------------------------------------------------

library(readr)     # fast CSV I/O
library(dplyr)     # data manipulation
library(tidyr)     # pivot_longer()
library(ggplot2)   # plotting
library(viridis)   # color scales (used inside plot helpers typically)
library(csthelpers) # scenario surface weight methods + plotting/eval helpers




# ------------------------------------------------------------------------------
# 1) Input data
# ------------------------------------------------------------------------------


# Lookup table defining the scenario grid and scenario metadata.
# Expected: columns include rlz (realization id), and at least tavg, prcp.
str_lookup <- read_csv("data/s4w_rhine_strlookup.csv")

# Use a single grid realization for plotting/computation.
# Assumption: rlz == 1 represents the baseline grid you want to evaluate on.
str_grid <- str_lookup  |>
  filter(rlz == 1)

# Ensemble summary stats for GCM samples (one row per model/scenario/horizon, etc.)
gcm_data <- read_csv("data/s4w_rhine_gcm_summary_stats.csv")

# Filter to a horizon and drop NAs (strict, but keeps downstream methods stable).
# Assumption: horizon == "far" exists and is the intended slice.
gcm_use <- gcm_data  |>
  filter(horizon == "far")  |>
  tidyr::drop_na()

# Impacts/output table for later weighting (used in the last section).
str_out <- read_csv("data/s4w_rhine_strtest_Q50.csv") |>
  pivot_longer(-strid, names_to = "location", values_to = "value")

# Plot limits derived from the grid (ensures comparable facets across methods).
xlim <- range(str_grid$tavg, na.rm = TRUE)
ylim <- range(str_grid$prcp, na.rm = TRUE)


# ------------------------------------------------------------------------------
# 3) Estimate scenario surface weights (3 methods)
# ------------------------------------------------------------------------------

# Common column mapping and method settings.
# Keep these centralized to avoid silent inconsistencies.
pr_col    <- "prcp"
ta_col    <- "tavg"
group_col <- "scenario"

# Common computational settings
area_weight <- "regular"  # "regular" typically means each grid cell equal area
scale       <- "global"   # global scaling of axes (per your package semantics)
normalize   <- TRUE       # normalize weights (e.g., sum to 1 per grid point)
verbose     <- TRUE

# --------------------------------------
# 3A) KDE (kernel density estimate)
# --------------------------------------

# KDE: nonparametric density estimate per scenario; flexible, can be data-hungry.
w_kde <- compute_scenario_surface_weights_kde(
  ensemble_data = gcm_use,
  scenario_grid = str_grid,
  pr_col        = pr_col,
  ta_col        = ta_col,
  group_col     = group_col,
  area_weight   = area_weight,
  scale         = scale,
  normalize     = normalize,
  verbose       = verbose
)

p_kde <- plot_kde(
  kde_grid   = w_kde,
  samples    = gcm_use,
  facet_col  = "scenario",
  z_col      = "weight",
  x_limits   = xlim,
  y_limits   = ylim
) + ggtitle("KDE")
p_kde

# --------------------------------------
# 3B) MVN (multivariate normal fit)
# --------------------------------------

# MVN: parametric (mean + covariance) per scenario; smooth and stable if fit is sane.
w_mvn <- compute_scenario_surface_weights_mvn(
  ensemble_data = gcm_use,
  scenario_grid = str_grid,
  pr_col        = pr_col,
  ta_col        = ta_col,
  group_col     = group_col,
  area_weight   = area_weight,
  scale         = scale,
  normalize     = normalize,
  verbose       = verbose
)

p_mvn  <- plot_kde(
  kde_grid   = w_mvn,
  samples    = gcm_use,
  facet_col  = "scenario",
  z_col      = "weight",
  x_limits   = xlim,
  y_limits   = ylim
) + ggtitle("MVN")
p_mvn

# --------------------------------------
# 3C) Gaussian copula
# --------------------------------------

# Copula: models dependence separately from marginals; can capture non-elliptical
# dependencies better than MVN, but Gaussian copula still implies Gaussian dependence
# structure in rank space (tail dependence limitations).
w_cop <- compute_scenario_surface_weights_cop(
  ensemble_data = gcm_use,
  scenario_grid = str_grid,
  pr_col        = pr_col,
  ta_col        = ta_col,
  group_col     = group_col,
  area_weight   = area_weight,
  scale         = scale,
  normalize     = normalize,
  verbose       = verbose
)


p_copula <- plot_kde(
  kde_grid   = w_cop,
  samples    = gcm_use,
  facet_col  = "scenario",
  z_col      = "weight",
  x_limits   = xlim,
  y_limits   = ylim
) + ggtitle("Copula")

p_copula

# ------------------------------------------------------------------------------
# 4) Evaluate scenario surface weights
# ------------------------------------------------------------------------------

# Evaluation uses the same mapping + scoring parameters for fairness.
# mapping = "knn" implies: at each GCM sample point, map it to the grid via KNN,
# then score how well the surface represents the ensemble (per your evaluator).

mapping       <- "knn"
k_neighbors   <- 5
tail_quantile <- 0.1  # lower/upper tails (10% and 90%) define "tail" regions

eval_kde <- evaluate_scenario_surface_weights(
  weight_surface = w_kde,
  ensemble_data  = gcm_use,
  ta_col         = ta_col,
  pr_col         = pr_col,
  group_col      = group_col,
  mapping        = mapping,
  k              = k_neighbors,
  tail_quantile  = tail_quantile,
  verbose        = TRUE
)

# BUGFIX: your original script overwrote eval_mvn and eval_cop by reassigning eval_kde.
# Keep them distinct.
eval_mvn <- evaluate_scenario_surface_weights(
  weight_surface = w_mvn,
  ensemble_data  = gcm_use,
  ta_col         = ta_col,
  pr_col         = pr_col,
  group_col      = group_col,
  mapping        = mapping,
  k              = k_neighbors,
  tail_quantile  = tail_quantile,
  verbose        = TRUE
)

eval_cop <- evaluate_scenario_surface_weights(
  weight_surface = w_cop,
  ensemble_data  = gcm_use,
  ta_col         = ta_col,
  pr_col         = pr_col,
  group_col      = group_col,
  mapping        = mapping,
  k              = k_neighbors,
  tail_quantile  = tail_quantile,
  verbose        = TRUE
)

# Compare the evaluation outputs side-by-side (assumes compare_methods() exists in csthelpers).
comparison <- compare_methods(
  kde    = eval_kde,
  mvn    = eval_mvn,
  copula = eval_cop
)

print(comparison)

# Raw comparison table (if compare_methods() returns a list with $comparison)
comparison$comparison

# ------------------------------------------------------------------------------
# 5) Apply weights to impacts for a site (example)
# ------------------------------------------------------------------------------

scn_wght <- str_lookup |>
  left_join(w_kde |> select(-strid, -rlz), by = c("tavg", "prcp"), relationship = "many-to-many")

# Custom columns
result <- compute_weighted_impacts(
  scenario_weights = scn_wght,
  impact_data = str_out_long,
  scenario_cols = list(scenario = c("scenario"), weight = "weight"),
  impact_cols = list(location = "location", value = "value"),
  realization_agg_fn = min, #stats::median,
  return_detail = TRUE
)

# Summary: median weighted impact per scenario (typically 4 rows for SSPs)
result$summary

# Detail: weighted impact per scenario per realization (dimensions depend on your inputs)
result$by_realization

