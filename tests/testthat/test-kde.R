# Functions: plot_kde (R/plot_kde.R)

testthat::test_that("plot_kde returns ggplot for valid inputs", {
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

  p <- plot_kde(kde_long, obs, bins = 3, raster_interpolate = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})
