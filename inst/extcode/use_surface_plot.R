

# Required libraries
library(readr)
library(dplyr)
library(ggplot2)

# Source script
source("R/climate_surface_plot.R")

# Climate stress test data
str_lookup <- read_csv("data/s4w_rhine_strlookup.csv")
str_out <- read_csv("data/s4w_rhine_strtest_Q50.csv") %>%
  left_join(str_lookup, by = "strid") %>%
  filter(rlz == 2) %>%
  select(strid, tavg, prcp, Lobith)

# GCM-Projections
gcm_data <- read_csv("data/s4w_rhine_gcm_summary_stats.csv")

data = str_out
x_var = "prcp"
y_var = "tavg"
z_var = "Lobith"
threshold = NULL
title = "Climate Response Surface"
x_label = expression(Delta ~ "Precipitation")
y_label = expression(Delta ~ "Temperature")
x_suffix = "%"
y_suffix = "\u00B0C"
failure_dir = 1
x_breaks = NULL
y_breaks = NULL
n_contours = 25
z_limits = NULL
panel_size_in = 6.0
legend_barheight_in = 0.25
text_size = 0.7
facet = FALSE
facet_by = NULL
facet_levels = NULL

range(str_out$Lobith)


# Report
legend_barwidth_spec <- 4.15
legend_barheight_spec <- 0.25
text_size <- 0.95
p_width <- 6.5
p_height <- 6


p1 <- climate_surface_base(
    data = str_out,
    x_var = x_var,
    y_var = y_var,
    z_var = z_var,
    threshold = threshold,
    title = title,
    x_label = x_label,
    y_label = y_label,
    x_suffix = x_suffix,
    y_suffix = y_suffix,
    failure_dir = failure_dir,
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    n_contours = n_contours,
    z_limits = z_limits,
    legend_barwidth_spec  = legend_barwidth_spec,
    legend_barheight_spec  = legend_barheight_spec,
    text_size = text_size,
    facet = facet,
    facet_by = facet_by,
    facet_levels = facet_levels)

ggplot2::ggsave("TEMP/figure1.png", p1, width = p_width, height = p_height, dpi = 300)


p2 <- climate_surface_gcm_overlay(
  p1, gcm_data = gcm_data,
  horizon_levels = "near",
  ellipse_levels = c(0.50, 0.90),
  ellipse_group = "none" # "none", "scenario", "horizon", "scenario_horizon"
)

ggplot2::ggsave("TEMP/figure2.png", p2, width = p_width, height = p_height, dpi = 300)

p3 <- climate_surface_gcm_overlay(
  p1, gcm_data = gcm_data,
  horizon_levels = "near",
  ellipse_levels = c(0.50, 0.90),
  ellipse_method = "robust",
  ellipse_group = "none"
)

ggplot2::ggsave("TEMP/figure3.png", p3, width = p_width, height = p_height, dpi = 300)


p4 <- climate_surface_gcm_overlay(
  p1, gcm_data = gcm_data,
  horizon_levels = "near",
  kde_levels = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
  kde_group = "none",
  kde_bw_adjust = 1.0,
  kde_n = 140

)
ggplot2::ggsave("TEMP/figure4.png", p4, width = p_width, height = p_height, dpi = 300)


