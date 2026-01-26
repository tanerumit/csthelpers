# tests/testthat/test-climate_surface_base.R

testthat::test_that("climate_surface_base() returns expected structure and ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("scales")

  # Simple grid with baseline point (0,0)
  x <- seq(-20, 20, by = 10)
  y <- seq(-2, 2, by = 1)
  df <- expand.grid(prcp = x, tavg = y)
  df$impact <- 100 + 2 * df$tavg - 0.5 * df$prcp
  # ensure baseline exists
  testthat::expect_true(any(df$prcp == 0 & df$tavg == 0))

  out <- climate_surface_base(
    data = df,
    x_var = "prcp",
    y_var = "tavg",
    z_var = "impact",
    title = "Test Surface"
  )

  testthat::expect_type(out, "list")
  testthat::expect_named(out, c("p", "width", "height"))
  testthat::expect_s3_class(out$p, "ggplot")
  testthat::expect_true(is.numeric(out$width) && length(out$width) == 1)
  testthat::expect_true(is.numeric(out$height) && length(out$height) == 1)

  # Metadata attached
  testthat::expect_true(is.numeric(attr(out$p, "threshold")))
  testthat::expect_identical(attr(out$p, "x_var"), "prcp")
  testthat::expect_identical(attr(out$p, "y_var"), "tavg")
  testthat::expect_identical(attr(out$p, "z_var"), "impact")
  testthat::expect_true(is.numeric(attr(out$p, "x_breaks")))
  testthat::expect_true(is.numeric(attr(out$p, "y_breaks")))
  testthat::expect_true(is.numeric(attr(out$p, "z_limits")))
  testthat::expect_true(is.numeric(attr(out$p, "contour_breaks")))
})

testthat::test_that("required arguments and column existence are validated", {
  df <- data.frame(a = 1:3, b = 1:3, c = 1:3)

  testthat::expect_error(
    climate_surface_base(data = NULL, x_var = "a", y_var = "b", z_var = "c"),
    "'data' must be a data.frame"
  )

  testthat::expect_error(
    climate_surface_base(data = df, x_var = NULL, y_var = "b", z_var = "c"),
    "'x_var', 'y_var', 'z_var' are required"
  )

  testthat::expect_error(
    climate_surface_base(data = df, x_var = "missing", y_var = "b", z_var = "c"),
    "Column 'missing' not found"
  )

  testthat::expect_error(
    climate_surface_base(data = df, x_var = "a", y_var = "missing", z_var = "c"),
    "Column 'missing' not found"
  )

  testthat::expect_error(
    climate_surface_base(data = df, x_var = "a", y_var = "b", z_var = "missing"),
    "Column 'missing' not found"
  )
})

testthat::test_that("x/y/z must be numeric", {
  df <- data.frame(
    prcp = c(0, 10, -10),
    tavg = c(0, 1, -1),
    impact = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  df_bad_x <- df; df_bad_x$prcp <- as.character(df_bad_x$prcp)
  testthat::expect_error(
    climate_surface_base(df_bad_x, "prcp", "tavg", "impact"),
    "'prcp' must be numeric"
  )

  df_bad_y <- df; df_bad_y$tavg <- as.character(df_bad_y$tavg)
  testthat::expect_error(
    climate_surface_base(df_bad_y, "prcp", "tavg", "impact"),
    "'tavg' must be numeric"
  )

  df_bad_z <- df; df_bad_z$impact <- as.character(df_bad_z$impact)
  testthat::expect_error(
    climate_surface_base(df_bad_z, "prcp", "tavg", "impact"),
    "'impact' must be numeric"
  )
})

testthat::test_that("failure_dir must be 1 or -1; n_contours must be >= 3", {
  df <- expand.grid(prcp = c(-10, 0, 10), tavg = c(-1, 0, 1))
  df$impact <- 1 + df$prcp + df$tavg

  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", failure_dir = 0),
    "'failure_dir' must be 1 or -1"
  )

  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", n_contours = 2),
    "n_contours.*>= 3"
  )
})

testthat::test_that("z_limits must be length-2, finite, and increasing; clipping is applied", {
  skip_if_not_installed("ggplot2")

  df <- expand.grid(prcp = seq(-20, 20, 10), tavg = seq(-2, 2, 1))
  df$impact <- 100 + 2 * df$tavg - 0.5 * df$prcp

  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", z_limits = c(1, 1)),
    "z_limits.*\\[1\\].*<.*\\[2\\]"
  )
  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", z_limits = c(NA, 2)),
    "z_limits"
  )

  out <- climate_surface_base(
    df, "prcp", "tavg", "impact",
    z_limits = c(90, 110),
    threshold = 100
  )

  # Metadata reports the legend/surface range
  testthat::expect_identical(attr(out$p, "z_limits"), c(90, 110))

  # The ggplot stores the (possibly clipped) data in p$data
  testthat::expect_true(all(out$p$data$impact >= 90 - 1e-12))
  testthat::expect_true(all(out$p$data$impact <= 110 + 1e-12))
})

testthat::test_that("threshold inference works when baseline exists, errors when it doesn't", {
  df <- expand.grid(prcp = c(-10, 0, 10), tavg = c(-1, 0, 1))
  df$impact <- 5 + df$prcp * 0 + df$tavg * 0
  # baseline exists and has impact 5
  out <- climate_surface_base(df, "prcp", "tavg", "impact")
  testthat::expect_equal(attr(out$p, "threshold"), 5)

  # remove baseline -> should error unless threshold provided
  df2 <- subset(df, !(prcp == 0 & tavg == 0))
  testthat::expect_error(
    climate_surface_base(df2, "prcp", "tavg", "impact", threshold = NULL),
    "Cannot infer threshold"
  )

  out2 <- climate_surface_base(df2, "prcp", "tavg", "impact", threshold = 5)
  testthat::expect_equal(attr(out2$p, "threshold"), 5)
})

testthat::test_that("z must have finite values and non-zero range", {
  df <- data.frame(prcp = c(-10, 0, 10), tavg = c(-1, 0, 1), impact = c(NA, NA, NA))
  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", threshold = 0),
    "z has no finite values"
  )

  df2 <- data.frame(prcp = c(-10, 0, 10), tavg = c(-1, 0, 1), impact = c(1, 1, 1))
  testthat::expect_error(
    climate_surface_base(df2, "prcp", "tavg", "impact", threshold = 1),
    "z has zero range"
  )
})

testthat::test_that("facet logic: requires facet_by; filters to facet_levels; adds facet_wrap", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- expand.grid(prcp = c(-10, 0, 10), tavg = c(-1, 0, 1), grp = c("A", "B"))
  df$impact <- 100 + 2 * df$tavg - 0.5 * df$prcp + ifelse(df$grp == "B", 1, 0)

  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", facet = TRUE, facet_by = NULL),
    "facet_by.*required"
  )
  testthat::expect_error(
    climate_surface_base(df, "prcp", "tavg", "impact", facet = TRUE, facet_by = "missing"),
    "Column 'missing' not found"
  )

  out <- climate_surface_base(
    df, "prcp", "tavg", "impact",
    facet = TRUE, facet_by = "grp", facet_levels = c("B")
  )

  # Data are filtered to facet_levels
  testthat::expect_true(all(out$p$data$grp %in% "B"))
  testthat::expect_s3_class(out$p$facet, "FacetWrap")
})

testthat::test_that("contour_breaks are within z_limits (or data range) and produce n_bins = length-1", {
  df <- expand.grid(prcp = seq(-20, 20, 10), tavg = seq(-2, 2, 1))
  df$impact <- 100 + 2 * df$tavg - 0.5 * df$prcp

  out <- climate_surface_base(
    df, "prcp", "tavg", "impact",
    threshold = 100,
    n_contours = 17,
    z_limits = c(90, 110)
  )

  br <- attr(out$p, "contour_breaks")
  testthat::expect_true(is.numeric(br))
  testthat::expect_true(all(is.finite(br)))
  testthat::expect_true(isTRUE(all(diff(br) > 0)))

  # Endpoints should be exactly z_limits (by design)
  testthat::expect_equal(br[1], 90)
  testthat::expect_equal(br[length(br)], 110)

  # n_bins used in filled contour is breaks-1; we at least ensure > 0
  testthat::expect_gt(length(br) - 1L, 0L)
})

testthat::test_that("x_breaks / y_breaks defaults are sorted unique of data, and axis limits respect them", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    prcp = c(10, 0, -10, 10),
    tavg = c(1, 0, -1, 1),
    impact = c(1, 2, 3, 4)
  )

  out <- climate_surface_base(df, "prcp", "tavg", "impact", threshold = 2.5)

  xb <- attr(out$p, "x_breaks")
  yb <- attr(out$p, "y_breaks")
  testthat::expect_identical(xb, sort(unique(df$prcp)))
  testthat::expect_identical(yb, sort(unique(df$tavg)))

  # Axis limits were set to range(x_breaks)/range(y_breaks) inside scales;
  # ggplot stores them in the scale objects. We only sanity-check that scales exist.
  sc_x <- out$p$scales$get_scales("x")
  sc_y <- out$p$scales$get_scales("y")
  testthat::expect_true(!is.null(sc_x))
  testthat::expect_true(!is.null(sc_y))
})
