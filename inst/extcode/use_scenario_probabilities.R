
# Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)

# Datsets
str_lookup <- read_csv("data/s4w_rhine_strlookup.csv")
str_grid <- str_lookup %>% filter(rlz == 1)

gcm_data <- read_csv("data/s4w_rhine_gcm_summary_stats.csv")
gcm_use <- gcm_data %>% filter(horizon == "far") %>% na.omit()
str_out <- read_csv("data/s4w_rhine_strtest_Q50.csv")
xlim <- range(str_grid$tavg)
ylim <- range(str_grid$prcp)

################################################################################

# KDE Based estimation
kde <- estimate_scenario_probs_kde(
  ensemble_data = gcm_use, scenario_grid = str_grid,
  pr_col = "prcp", ta_col = "tavg", group_col = "scenario",
  area_weight = "regular", scale = "global", normalize = TRUE, verbose = TRUE)

p <- plot_kde_weights(kde = kde, obs_data = gcm_use,
        xlim = xlim, ylim = ylim); p


weighted_lobith <- compute_weighted_impacts_site(
  scenario_key = str_lookup,   # or your scenario_key tibble
  kde_weights = kde,
  impacts      = str_out,
  site_col     = "Lobith",
  return_detail = TRUE
)

weighted_lobith$summary        # 4 rows: median weighted impact per scenario
weighted_lobith$by_realization # 4 x 5 rows: weighted impact per scenario per rlz




