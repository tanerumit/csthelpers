

library(ggplot2)
library(dplyr)
library(sf)


# S4W - Rhine Basin example

rivers_sf <- read_sf("C:/Users/taner/WS/s4w/data/staticgeoms/rivers.geojson")
basins_sf <- read_sf("C:/Users/taner/WS/s4w/data/staticgeoms/basins.geojson")
subbasins_sf <- read_sf("C:/Users/taner/WS/s4w/data/staticgeoms/subcatch_wflow-gauges-ahr.geojson")
gauges_sf <- read_sf("C:/Users/taner/WS/s4w/data/staticgeoms/gauges_grade_climate.geojson")

p <- plot_basin_map(
  basins_sf = basins_sf,
  rivers_sf = rivers_sf,
  subbasins_sf = subbasins_sf,
  points_sf = gauges_sf,
  label_col = "name",
  stream_order_col = "strord",
  min_stream_order = 3,
  bbox_buffer = 0.5,
  add_scale = TRUE,
  add_north = TRUE
)

p
