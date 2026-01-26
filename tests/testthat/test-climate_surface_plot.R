# Functions: climate_surface_base (R/climate_surface_plot.R), ggsave_smart (R/climate_surface_plot.R)

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
  testthat::expect_true(is.numeric(attr(p, "legend_barwidth_in")))
  testthat::expect_true(is.numeric(attr(p, "legend_barheight_in")))
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

testthat::test_that("ggsave_smart writes output and returns sizing details", {
  testthat::skip_if_not_installed("ggplot2")
  suppressPackageStartupMessages(library(ggplot2))

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()

  out_file <- tempfile(fileext = ".png")
  res <- ggsave_smart(out_file, p, target = "ppt", verbose = FALSE)

  testthat::expect_true(file.exists(out_file))
  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("width_in", "height_in", "dpi") %in% names(res)))
})
