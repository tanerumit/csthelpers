
library(csthelpers)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# Stress testing grid (base surface)
str_data <- readr::read_csv("data/Qstats.csv") %>%
  filter(statistic == "mean") %>%
  select(tavg, prcp, value = Q_1)

p_base <- climate_surface_base(
  data = str_data,
  x_var = "prcp",
  y_var = "tavg",
  z_var = "value",
  threshold = NULL,
  title = "Response Surface",
  x_label = "Changes in Annual Precipitation" ,
  y_label = "Changes in Annual Temperature" ,
  x_suffix = "%",
  y_suffix = "\u00B0C",
  failure_dir = 1,
  x_breaks = NULL,
  y_breaks = NULL,
  n_contours = 11,
  z_limits = NULL, ##c(0,1300),
  panel_size_in = 6.0,
  legend_barheight_in = 0.25,
  text_size = 0.7,
  # Faceting
  facet = FALSE,
  facet_by = NULL,
  facet_levels = NULL); p_base


ggsave("examples/surface.png", p_base$p, width = p_base$width, height = p_base$height,
       units = "in", dpi = 300)

gcm_df <- readr::read_csv("data/annual_change_scalar_stats_summary_mean.csv") %>%
  dplyr::filter(stats == "mean")  #

p_gcm <- climate_surface_gcm_overlay(
  p = p_base$p,
  gcm_data = gcm_df,
  x_var = "prcp",
  y_var = "tavg",
  scenario_levels = NULL,
  horizon_levels  = "far",
  ellipse_levels = c(0.50, 0.90),   # optional
  ellipse_group = "scenario")

ggsave("examples/surface2.png",
       plot = p_gcm$p, width = p_gcm$width, height =  p_gcm$height,
        units = "in", dpi = 300)




################################################################################


# Stress testing grid (base surface)
str_data <- readr::read_csv("data/Qstats.csv") %>%
  filter(statistic == "mean") %>%
  pivot_longer(cols = Q_1:Q_1622500, names_to = "station", values_to = "value")


p3 <- climate_surface_base(
  data = str_data,
  x_var = "prcp",
  y_var = "tavg",
  z_var = "value",
  threshold = NULL,
  title = "Response Surface",
  x_label = "Changes in Annual Precipitation" ,
  y_label = "Changes in Annual Temperature" ,
  x_suffix = "%",
  y_suffix = "\u00B0C",
  failure_dir = 1,
  x_breaks = NULL,
  y_breaks = NULL,
  n_contours = 11,
  z_limits = NULL,
  panel_size_in = 6.0,
  legend_barheight_in = 0.25,
  text_size = 0.7,
  # Faceting
  facet = TRUE,
  facet_by = "station",
  facet_levels = c("Q_1", "Q_4003", "Q_4005", "Q_5001"))

ggsave("examples/surface2.png", plot = p3$p, width = 8, height =  8.5)
