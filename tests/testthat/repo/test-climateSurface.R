
# Functions: climate_surface_base (R/plot_climate_surface.R),
#   climate_surface_gcm_overlay (R/plot_climate_surface.R)

test_that("climate_surface_base returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("scales")

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 5),
    temp_change = seq(-2, 4, by = 2)
  )
  df$response <- with(df, 0.5 + 0.01 * prcp_change - 0.02 * temp_change)

  p <- climate_surface_base(
    data = df,
    x_var = "prcp_change",
    y_var = "temp_change",
    z_var = "response",
    threshold = 0.5
  )

  expect_s3_class(p, "ggplot")
  expect_true(is.list(p$layers))
  expect_true(length(p$layers) >= 1)
  expect_identical(attr(p, "x_var"), "prcp_change")
  expect_identical(attr(p, "y_var"), "temp_change")
})

test_that("threshold is used and produces a contour layer", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("scales")

  df <- expand.grid(
    prcp_change = seq(-5, 5, by = 1),
    temp_change = seq(-1, 3, by = 1)
  )
  df$response <- with(df, 0.8 + 0.02 * prcp_change + 0.05 * temp_change)

  p <- climate_surface_base(
    data = df,
    x_var = "prcp_change",
    y_var = "temp_change",
    z_var = "response",
    threshold = 0.9
  )

  geom_classes <- vapply(p$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomContour" %in% geom_classes)
})

test_that("custom breaks are stored on the plot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("scales")

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 5),
    temp_change = seq(-2, 4, by = 2)
  )
  df$response <- runif(nrow(df))

  x_breaks <- c(-10, 0, 10)
  y_breaks <- c(-2, 1, 4)

  p <- climate_surface_base(
    data = df,
    x_var = "prcp_change",
    y_var = "temp_change",
    z_var = "response",
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    threshold = 0.5
  )

  expect_identical(attr(p, "x_breaks"), x_breaks)
  expect_identical(attr(p, "y_breaks"), y_breaks)
})

test_that("GCM overlay does not error and adds points", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("scales")

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 5),
    temp_change = seq(-2, 6, by = 2)
  )
  df$response <- runif(nrow(df))

  gcm <- data.frame(
    prcp_change = c(-5, 0, 5),
    temp_change = c(0, 2, 4),
    scenario = c("ssp245", "ssp585", "ssp245"),
    horizon  = c("near", "far", "near")
  )

  p_base <- climate_surface_base(
    data = df,
    x_var = "prcp_change",
    y_var = "temp_change",
    z_var = "response",
    threshold = 0.5
  )

  p <- climate_surface_gcm_overlay(
    p = p_base,
    gcm_data = gcm,
    scenario_levels = c("ssp245", "ssp585"),
    horizon_levels = c("near", "far"),
    spread_method = "none",
    show_legend = FALSE
  )

  geom_classes <- vapply(p$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomPoint" %in% geom_classes)
})

test_that("facet filtering works and returns facets", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("scales")

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 10),
    temp_change = seq(-2, 6, by = 2),
    basin = c("A", "B")
  )
  df$response <- runif(nrow(df))

  p <- climate_surface_base(
    data = df,
    x_var = "prcp_change",
    y_var = "temp_change",
    z_var = "response",
    threshold = 0.5,
    facet = TRUE,
    facet_by = "basin",
    facet_levels = c("A", "B")
  )

  expect_s3_class(p$facet, "FacetWrap")
  expect_true(all(p$data$basin %in% c("A", "B")))
})
