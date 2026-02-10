test_that("radial_plot() rejects non-data.frame or empty data", {
  expect_error(radial_plot(NULL, category = a, value = b), "`data` must be a non-empty data.frame")
  expect_error(radial_plot(list(a = 1), category = a, value = b), "`data` must be a non-empty data.frame")
  expect_error(radial_plot(data.frame(), category = a, value = b), "`data` must be a non-empty data.frame")
})

test_that("radial_plot() validates scalar controls", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))

  expect_error(radial_plot(df, cat, val, start_angle = "90"), "`start_angle` must be a single non-NA numeric")
  expect_error(radial_plot(df, cat, val, start_angle = c(90, 91)), "`start_angle` must be a single non-NA numeric")

  expect_error(radial_plot(df, cat, val, radial_expand = 0.99), "`radial_expand` must be between 1 and 2")
  expect_error(radial_plot(df, cat, val, radial_expand = 2.01), "`radial_expand` must be between 1 and 2")

  expect_error(radial_plot(df, cat, val, clockwise = 1), "`clockwise` must be a single non-NA logical")
  expect_error(radial_plot(df, cat, val, show_spokes = "yes"), "`show_spokes` must be a single non-NA logical")
  expect_error(radial_plot(df, cat, val, show_grid_circles = NA), "`show_grid_circles` must be a single non-NA logical")
  expect_error(radial_plot(df, cat, val, show_value_labels = 0), "`show_value_labels` must be a single non-NA logical")
})

test_that("radial_plot() errors when category/value columns are missing", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))

  expect_error(radial_plot(df, missing_cat, val))
  expect_error(radial_plot(df, cat, missing_val))
})

test_that("radial_plot() enforces numeric value column", {
  df <- data.frame(cat = c("A", "B"), val = c("1", "2"))
  expect_error(radial_plot(df, cat, val), "`value` column must be numeric")
})

test_that("radial_plot() handles NA values per na_handling", {
  df <- data.frame(cat = c("A", "B", "C"), val = c(1, NA_real_, 3))

  expect_warning(radial_plot(df, cat, val, na_handling = "warn"))
  expect_silent(radial_plot(df, cat, val, na_handling = "silent"))
  expect_error(radial_plot(df, cat, val, na_handling = "error"), "`value` contains 1 NA\\(s\\)")

  df_all_na <- data.frame(cat = c("A", "B"), val = c(NA_real_, NA_real_))
  expect_error(radial_plot(df_all_na, cat, val), "`value` contains only NA values")
})

test_that("radial_plot() validates radial_breaks", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))

  expect_error(radial_plot(df, cat, val, radial_breaks = "x"), "`radial_breaks` must be numeric")
  expect_error(radial_plot(df, cat, val, radial_breaks = c(0, NA)), "`radial_breaks` must be numeric, length >= 2, no NAs")
  expect_error(radial_plot(df, cat, val, radial_breaks = c(1, 2)), "`radial_breaks` must include 0")

  p <- radial_plot(df, cat, val, radial_breaks = c(0, 1, 2))
  expect_s3_class(p, "ggplot")
})

test_that("radial_plot() works with facet as NULL and as a column", {
  df <- data.frame(
    cat = rep(c("A", "B"), 2),
    val = c(1, 2, 3, 4),
    scen = rep(c("S1", "S2"), each = 2)
  )

  p1 <- radial_plot(df, cat, val, facet = NULL)
  expect_s3_class(p1, "ggplot")

  p2 <- radial_plot(df, cat, val, facet = scen)
  expect_s3_class(p2, "ggplot")
})

test_that("radial_plot() accepts column names passed as strings", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))
  p <- radial_plot(df, category = "cat", value = "val")
  expect_s3_class(p, "ggplot")
})

test_that("fill_mode = none returns ggplot", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))
  p <- radial_plot(df, cat, val, fill_mode = "none")
  expect_s3_class(p, "ggplot")
})

test_that("fill_mode = binned validates fill_levels and fill_colors", {
  df <- data.frame(cat = c("A", "B", "C"), val = c(1, 5, 9))

  expect_error(
    radial_plot(df, cat, val, fill_mode = "binned", fill_levels = NULL),
    "For `fill_mode = \"binned\"`, `fill_levels` must be a numeric vector"
  )
  expect_error(
    radial_plot(df, cat, val, fill_mode = "binned", fill_levels = c(0, NA)),
    "For `fill_mode = \"binned\"`, `fill_levels` must be a numeric vector"
  )
  expect_error(
    radial_plot(df, cat, val, fill_mode = "binned", fill_levels = c(1, 1, 1)),
    "`fill_levels` must contain at least 2 unique breakpoints"
  )

  # unnamed fill_colors length mismatch
  expect_error(
    radial_plot(df, cat, val, fill_mode = "binned", fill_levels = c(0, 5, 10), fill_colors = c("red")),
    "Unnamed `fill_colors` must have length 2"
  )

  # named fill_colors length must match bins
  expect_error(
    radial_plot(
      df, cat, val,
      fill_mode = "binned",
      fill_levels = c(0, 5, 10),
      fill_colors = c(Low = "red")
    ),
    "For binned fill, if `fill_colors` is named it must have exactly one color per bin"
  )

  p <- radial_plot(df, cat, val, fill_mode = "binned", fill_levels = c(0, 5, 10))
  expect_s3_class(p, "ggplot")

  # named colors aligned via names as labels (preferred path)
  p2 <- radial_plot(
    df, cat, val,
    fill_mode = "binned",
    fill_levels = c(0, 5, 10),
    fill_colors = c(Low = "green", High = "red")
  )
  expect_s3_class(p2, "ggplot")
})

test_that("fill_mode = continuous validates fill_colors and fill_levels", {
  df <- data.frame(cat = c("A", "B", "C"), val = c(1, 5, 9))

  expect_error(
    radial_plot(df, cat, val, fill_mode = "continuous", fill_colors = "red"),
    "For `fill_mode = \"continuous\"`, `fill_colors` must be a character vector with length >= 2"
  )

  expect_error(
    radial_plot(df, cat, val, fill_mode = "continuous", fill_levels = c(0, NA)),
    "For continuous fill, `fill_levels` must be numeric length >= 2 with no NAs"
  )
  expect_error(
    radial_plot(df, cat, val, fill_mode = "continuous", fill_levels = c(0, Inf)),
    "For continuous fill, `fill_levels` must be finite"
  )
  expect_error(
    radial_plot(df, cat, val, fill_mode = "continuous", fill_levels = c(1, 1)),
    "must contain at least 2 unique stop points"
  )

  # Works with default colors
  p <- radial_plot(df, cat, val, fill_mode = "continuous")
  expect_s3_class(p, "ggplot")

  # Works with >2 stops (forces interpolation)
  p2 <- radial_plot(
    df, cat, val,
    fill_mode = "continuous",
    fill_colors = c("yellow", "red"),
    fill_levels = c(0, 3, 6, 9)
  )
  expect_s3_class(p2, "ggplot")
})

test_that("invalid choices are rejected by match.arg", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))

  expect_error(radial_plot(df, cat, val, spokes_extend_to = "bad"))
  expect_error(radial_plot(df, cat, val, theme_style = "bad"))
  expect_error(radial_plot(df, cat, val, na_handling = "bad"))
  expect_error(radial_plot(df, cat, val, fill_mode = "bad"))
})

test_that("style argument: non-radial_style is replaced with default", {
  df <- data.frame(cat = c("A", "B"), val = c(1, 2))
  p <- radial_plot(df, cat, val, style = list(foo = 1))
  expect_s3_class(p, "ggplot")
})


