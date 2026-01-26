
library(testthat)
library(ggplot2)

test_that("climate_surface returns a ggplot object", {

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 5),
    temp_change = seq(-2, 4, by = 2)
  )
  df$response <- with(df, 0.5 + 0.01 * prcp_change - 0.02 * temp_change)

  p <- climate_surface(
    str.data   = df,
    variable.x = "prcp_change",
    variable.y = "temp_change",
    variable.z = "response",
    threshold.z = 0.5
  )

  expect_s3_class(p, "ggplot")
  expect_true(is.list(p$layers))
  expect_true(length(p$layers) >= 1)
})

test_that("threshold.z is used and produces a contour layer", {

  df <- expand.grid(
    prcp_change = seq(-5, 5, by = 1),
    temp_change = seq(-1, 3, by = 1)
  )
  df$response <- with(df, 0.8 + 0.02 * prcp_change + 0.05 * temp_change)

  p <- climate_surface(
    str.data   = df,
    variable.x = "prcp_change",
    variable.y = "temp_change",
    variable.z = "response",
    threshold.z = 0.9
  )

  geom_classes <- vapply(p$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomContour" %in% geom_classes)
})

test_that("custom breaks are used when provided", {

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 5),
    temp_change = seq(-2, 4, by = 2)
  )
  df$response <- runif(nrow(df))

  x_breaks <- c(-10, 0, 10)
  y_breaks <- c(-2, 1, 4)

  p <- climate_surface(
    str.data = df,
    variable.x = "prcp_change",
    variable.y = "temp_change",
    variable.z = "response",
    variable.x.breaks = x_breaks,
    variable.y.breaks = y_breaks
  )

  x_scale <- p$scales$get_scales("x")
  y_scale <- p$scales$get_scales("y")

  expect_identical(x_scale$breaks, x_breaks)
  expect_identical(y_scale$breaks, y_breaks)
})

test_that("GCM overlay does not error and adds points", {

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

  p <- climate_surface(
    str.data = df,
    gcm.data = gcm,
    variable.x = "prcp_change",
    variable.y = "temp_change",
    variable.z = "response"
  )

  geom_classes <- vapply(p$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomPoint" %in% geom_classes)
})

test_that("multi.panel filtering works and returns facets", {

  df <- expand.grid(
    prcp_change = seq(-10, 10, by = 10),
    temp_change = seq(-2, 6, by = 2),
    basin = c("A", "B")
  )
  df$response <- runif(nrow(df))

  p <- climate_surface(
    str.data = df,
    variable.x = "prcp_change",
    variable.y = "temp_change",
    variable.z = "response",
    multi.panel = TRUE,
    panel.variable = "basin",
    panel.variable.levels = c("A", "B")
  )

  expect_true("facet_wrap" %in% class(p$facet))
  expect_identical(p$facet$params$facets, ~ basin)
})
