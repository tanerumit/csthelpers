

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


##################### DATA/INPUT

# GCM projections for each SSP
gcm.data <- readr::read_csv("data/annual_change_scalar_stats_summary_mean.csv") %>%
  filter(horizon == "near")

# Climate stress testing grid
str.data <- readr::read_csv("data/Qstats.csv") %>%
  filter(statistic == "mean") %>% select(tavg, prcp, value = Q_1)
#readr::write_csv(stress_grid, "stress_test_grid.csv")


p <- climate_surface(
  str.data = str_data,
  gcm.data = gcm_data,
  variable.x = "tavg",
  variable.y = "prcp",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = "change in variable",
  variable.x.label = expression(Delta ~ "Precipitation"),
  variable.y.label = expression(Delta ~ "Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = FALSE,
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.z.breaks = NULL,
  variable.x.breaks = NULL,
  variable.y.breaks = NULL,
  contour.num = 15,
  text.scale = 0.6,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)
p
