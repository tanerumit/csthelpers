# Functions: climate_surface_base (R/plot_climate_surface.R), climate_surface_gcm_overlay (R/plot_climate_surface.R)

if (!exists("climate_surface_base", mode = "function") ||
    !exists("climate_surface_gcm_overlay", mode = "function")) {
  source_candidates <- c(
    file.path("R", "plot_climate_surface.R"),
    file.path("..", "..", "R", "plot_climate_surface.R")
  )
  source_file <- source_candidates[file.exists(source_candidates)][1]
  if (!is.na(source_file)) {
    sys.source(source_file, envir = environment())
  }
}

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
  testthat::expect_true(is.numeric(attr(p, "threshold")))
  testthat::expect_true(is.numeric(attr(p, "contour_breaks")))
  testthat::expect_true(all(diff(attr(p, "contour_breaks")) > 0))
  testthat::expect_equal(
    attr(p, "contour_breaks"),
    sort(unique(attr(p, "contour_breaks")))
  )
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
  geom_classes <- vapply(p_overlay$layers, function(x) class(x$geom)[1], character(1))
  testthat::expect_true("GeomPoint" %in% geom_classes)
})

testthat::test_that("surface_style overrides baseline_tol and legend sizing", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("scales")
  suppressPackageStartupMessages(library(ggplot2))

  df <- data.frame(
    dx = c(1e-5, 1e-5, 1, 1),
    dy = c(1e-5, 1, 1e-5, 1),
    dz = c(0.2, 0.3, 0.4, 0.5)
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

  sty <- surface_style(
    baseline_tol = 1e-4,
    legend_barwidth_spec = 5.5,
    legend_barheight_spec = 0.4,
    legend_max_labels = 2L
  )

  p <- climate_surface_base(
    data = df,
    x_var = "dx",
    y_var = "dy",
    z_var = "dz",
    threshold = NULL,
    n_contours = 5,
    style = sty
  )

  testthat::expect_equal(attr(p, "threshold"), 0.2)
  testthat::expect_equal(attr(p, "legend_barwidth_spec"), 5.5)
  testthat::expect_equal(attr(p, "legend_barheight_spec"), 0.4)
})

testthat::test_that("surface_style validates legend_max_labels >= 2", {
  testthat::expect_error(
    surface_style(legend_max_labels = 1L),
    "legend_max_labels.*>= 2"
  )
})

testthat::test_that("climate_surface_gcm_overlay validates alpha, size, and stroke", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("scales")
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

  testthat::expect_error(
    climate_surface_gcm_overlay(
      p = p_base,
      gcm_data = gcm_data,
      alpha = 1.5
    ),
    "alpha.*between 0 and 1"
  )

  testthat::expect_error(
    climate_surface_gcm_overlay(
      p = p_base,
      gcm_data = gcm_data,
      size = -0.1
    ),
    "size.*between 0 and Inf"
  )

  testthat::expect_error(
    climate_surface_gcm_overlay(
      p = p_base,
      gcm_data = gcm_data,
      stroke = -0.1
    ),
    "stroke.*between 0 and Inf"
  )
})
