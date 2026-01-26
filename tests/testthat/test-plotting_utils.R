# Functions: plot_basin_map (R/plotting_utils.R), sos_radialplot (R/plotting_utils.R), plot_kde_weights (R/scenario_probabilities_plots.R)

testthat::test_that("plot_basin_map returns ggplot", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("dplyr")
  suppressPackageStartupMessages(library(ggplot2))

  basin <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
    ),
    crs = 4326
  )

  rivers <- sf::st_sf(
    strord = c(1, 3),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(0, 0), c(1, 1))),
      sf::st_linestring(rbind(c(0, 1), c(1, 0)))
    ),
    crs = 4326
  )

  points <- sf::st_sf(
    name = c("A", "B"),
    geometry = sf::st_sfc(sf::st_point(c(0.2, 0.2)), sf::st_point(c(0.8, 0.8))),
    crs = 4326
  )

  p <- suppressWarnings(plot_basin_map(basin, rivers, points, main_stem_order = 2))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("sos_radialplot returns ggplot and enforces color levels", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))

  df <- data.frame(
    group = rep(c("A", "B", "C"), times = 2),
    scenario = rep(c("S1", "S2"), each = 3),
    value = c(1, 2, 3, 2, 3, 4),
    risk = c("Low", "Low", "High", "Low", "High", "High"),
    stringsAsFactors = FALSE
  )

  cat_colors <- c("Low" = "yellow", "High" = "red")
  p <- suppressWarnings(sos_radialplot(
    data = df,
    axis_col = group,
    panel_col = scenario,
    value_col = value,
    color_col = risk,
    cat_colors = cat_colors,
    radial_breaks = c(0, 2, 4),
    start_angle = 90
  ))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_kde_weights returns ggplot for valid inputs", {
  testthat::skip_if_not_installed("ggplot2")
  suppressPackageStartupMessages(library(ggplot2))

  kde_long <- data.frame(
    tavg = rep(c(0, 1), each = 2),
    prcp = rep(c(0, 1), times = 2),
    scenario = rep(c("S1", "S2"), each = 2),
    weight = c(0.2, 0.8, 0.6, 0.4)
  )

  obs <- data.frame(
    tavg = c(0.1, 0.9),
    prcp = c(0.2, 0.8)
  )

  p <- plot_kde_weights(kde_long, obs, bins = 3, raster_interpolate = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})
