
# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


##################### DATA/INPUT

# GCM projections for each SSP
gcm_data <- readr::read_csv("data/annual_change_scalar_stats_summary_mean.csv")

# Climate stress testing grid
str_data <- readr::read_csv("data/Qstats.csv") %>% filter(statistic == "mean")
stress_grid <- str_data %>% select(tavg, prcp, value = Q_1)
#readr::write_csv(stress_grid, "stress_test_grid.csv")

# Genealogy data from Kuma et al paper
kuma <- readr::read_csv("data/2022ms003588.csv")


################### CALCULATIONS
gcm_df <- gcm_data %>% filter(horizon == "far")

gcm_dfw <- compute_weights_genealogy(
  gcm_data   = gcm_df,
  kuma_table = kuma,
  model_col = "model",
  scenario_col = "scenario",
  cmip_phase = "CMIP6",
  method = "family",
  verbose = TRUE)

# Estimate Kernel Densities
prob_kde_fine <- estimate_scenario_probabilities_kde(
  ensemble_data = gcm_dfw,
  scenario_grid = stress_grid,
  weight_col = "w_genealogy",
  pr_col = "prcp",
  ta_col = "tavg",
  group_col = "scenario",
  bw = NULL,
  k = c(1.5, 2.0),
  alpha = 1.0,
  bw_min = c(0, 0),
  bw_max = c(Inf, Inf),
  min_samples = 5L,
  normalize = TRUE,
  chunk_size = 5000L,
  verbose = TRUE)

df_plot <- prob_kde_fine

p <- ggplot(df_plot) +
  theme_bw() +
  geom_raster(aes(tavg, prcp, fill = ssp585), interpolate = TRUE) +
  geom_contour(aes(tavg, prcp, z = ssp585),
    bins = 15, color = "black", linejoin = "round", lineend  = "round", alpha = 0.9,
    linewidth = 0.5) +
  scale_fill_viridis_c(
    name = "Relative\nweight",
    guide = guide_colorbar(barheight = unit(200, "pt"), ticks = TRUE)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

p

# 1) Estimate probabilities over stress-test grid
grid_probs <- estimate_grid_probabilities(
  gcm_df_w = gcm_dfw,
  stress_test_grid = stress_grid,
  scenario_col = "scenario",
  weight_col = "w_genealogy",
  bw = NULL    # auto bandwidth
)

# 2) Report weighted impact statistics
thresholds <- list(
  impact_gt_20 = list(value = 20, direction = "gt")
)

impact_stats <- report_weighted_impacts_from_grid(
  grid_probs = grid_probs,
  stress_test_grid = stress_grid,
  value_col = "value",
  thresholds = thresholds
)


