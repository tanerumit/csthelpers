# Functions: climate_surface_base (R/plot_climate_surface.R), climate_surface_gcm_overlay (R/plot_climate_surface.R)

testthat::test_that("climate_surface_base returns ggplot with metadata", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("scales")
  suppressPackageStartupMessages(library(ggplot2))

  df <- data.frame(
    dx = c(0, 0, 1, 1),
    dy = c(0, 1, 0, 1),
    dz = c(0.1, 0.2, 0.3, 0.4)
  )

  p <- climate_surface_base(
    data = df,
    x_var = "dx",
    y_var = "dy",
    z_var = "dz",
    threshold = 0.2,
    n_contours = 5
  )

  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(is.numeric(attr(p, "legend_barwidth_spec")))
  testthat::expect_true(is.numeric(attr(p, "legend_barheight_spec")))
})

testthat::test_that("climate_surface_base errors when baseline missing and threshold NULL", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("scales")
  suppressPackageStartupMessages(library(ggplot2))

  df <- data.frame(
    dx = c(1, 1),
    dy = c(1, 2),
    dz = c(0.1, 0.2)
  )

  testthat::expect_error(
    climate_surface_base(
      data = df,
      x_var = "dx",
      y_var = "dy",
      z_var = "dz",
      threshold = NULL,
      n_contours = 5
    ),
    "Cannot infer threshold"
  )
})

testthat::test_that("climate_surface_gcm_overlay adds scatter layer", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  suppressPackageStartupMessages(library(ggplot2))

  df <- data.frame(
    dx = c(0, 0, 1, 1),
    dy = c(0, 1, 0, 1),
    dz = c(0.1, 0.2, 0.3, 0.4)
  )

  p_base <- climate_surface_base(
    data = df,
    x_var = "dx",
    y_var = "dy",
    z_var = "dz",
    threshold = 0.2,
    n_contours = 5
  )

  gcm_data <- data.frame(
    scenario = c("ssp126", "ssp245"),
    horizon = c("near", "far"),
    dx = c(0.2, 0.8),
    dy = c(0.3, 0.7)
  )

  p_overlay <- climate_surface_gcm_overlay(
    p = p_base,
    gcm_data = gcm_data,
    scenario_levels = c("ssp126", "ssp245"),
    horizon_levels = c("near", "far"),
    spread_method = "none",
    show_legend = FALSE
  )

  testthat::expect_s3_class(p_overlay, "ggplot")
})
