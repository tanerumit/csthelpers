


# Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(csthelpers)

# Datasets
str_lookup <- read_csv("data/s4w_rhine_strlookup.csv")
str_grid <- str_lookup %>% filter(rlz == 1)

gcm_data <- read_csv("data/s4w_rhine_gcm_summary_stats.csv")
gcm_use <- gcm_data %>% filter(horizon == "far") %>% na.omit()
str_out <- read_csv("data/s4w_rhine_strtest_Q50.csv")
xlim <- range(str_grid$tavg)
ylim <- range(str_grid$prcp)

################################################################################

# Kernel Density based estimation
w_kde <- compute_scenario_surface_weights_kde(
  ensemble_data = gcm_use,
  scenario_grid = str_grid,
  pr_col = "prcp",
  ta_col = "tavg",
  group_col = "scenario",
  area_weight = "regular",
  scale = "global",
  normalize = TRUE,
  verbose = TRUE)

w_kde_long <- w_kde %>%
  pivot_longer(ssp126:ssp585, names_to = "scenario", values_to = "value")

p_kde <- plot_kde(kde_grid = w_kde_long, samples = gcm_use,
              facet_col = "scenario",
              z_col = "value",
              x_limits = xlim, y_limits = ylim) +
              ggtitle("KDE"); p_kde

# Estimate grid weights using multivariate normal
w_mvn <- compute_scenario_surface_weights_mvnorm(
  ensemble_data = gcm_use,
  scenario_grid = str_grid,
  pr_col = "prcp",
  ta_col = "tavg",
  group_col = "scenario",
  area_weight = "regular",
  scale = "global",
  normalize = TRUE,
  verbose = TRUE)


w_mvn_long <- w_mvn %>%
  pivot_longer(ssp126:ssp585, names_to = "scenario", values_to = "value")


p_mvnorm <- plot_kde(kde_grid = w_mvn_long, samples = gcm_use,
              facet_col = "scenario",
              z_col = "value",
              x_limits = xlim, y_limits = ylim) +
              ggtitle("MVNORM"); p_mvnorm



# Estimate using gaussian copula
w_cop <- compute_scenario_surface_weights_copula(
  ensemble_data = gcm_use,
  scenario_grid = str_grid,
  pr_col = "prcp",
  ta_col = "tavg",
  group_col = "scenario",
  area_weight = "regular",
  scale = "global",
  normalize = TRUE,
  verbose = TRUE)

w_cop_long <- w_cop %>%
  pivot_longer(ssp126:ssp585, names_to = "scenario", values_to = "value")


p_copula <- plot_kde(kde_grid = w_cop_long, samples = gcm_use,
                     facet_col = "scenario",
                     z_col = "value",
                     x_limits = xlim, y_limits = ylim) +
                     ggtitle("Copula"); p_copula


#################################################################################

# score surfaces (same evaluator)
eval_kde <- evaluate_scenario_surface_weights(
  weight_surface = w_kde,
  ensemble_data = gcm_use,
  ta_col = "tavg",
  pr_col = "prcp",
  group_col = "scenario",
  mapping = "knn",
  k = 5,
  tail_quantile = 0.1,  # 10th/90th percentile defines "tail"
  verbose = TRUE
)

eval_mvn <- eval_kde <- evaluate_scenario_surface_weights(
  weight_surface = w_mvn,
  ensemble_data = gcm_use,
  ta_col = "tavg",
  pr_col = "prcp",
  group_col = "scenario",
  mapping = "knn",
  k = 5,
  tail_quantile = 0.1,  # 10th/90th percentile defines "tail"
  verbose = TRUE
)

eval_cop <- eval_kde <- evaluate_scenario_surface_weights(
  weight_surface = w_cop,
  ensemble_data = gcm_use,
  ta_col = "tavg",
  pr_col = "prcp",
  group_col = "scenario",
  mapping = "knn",
  k = 5,
  tail_quantile = 0.1,  # 10th/90th percentile defines "tail"
  verbose = TRUE
)

comparison <- compare_methods(
  kde = eval_kde,
  mvn = eval_mvn,
  copula = eval_cop
)

print(comparison)

# Access raw comparison table
comparison$comparison


################################################################################
################################################################################
################################################################################
################################################################################

weighted_lobith <- compute_weighted_impacts_site(
  scenario_key = str_lookup,   # or your scenario_key tibble
  kde_weights = kde,
  impacts      = str_out,
  site_col     = "Lobith",
  return_detail = TRUE
)

weighted_lobith$summary        # 4 rows: median weighted impact per scenario
weighted_lobith$by_realization # 4 x 5 rows: weighted impact per scenario per rlz




