layer_geom_classes <- function(plot) {
  vapply(plot$layers, function(layer) class(layer$geom)[1], character(1))
}

fill_scale <- function(plot) {
  plot$scales$get_scales("fill")
}

y_scale <- function(plot) {
  plot$scales$get_scales("y")
}

test_that("radial_style() returns the expected structure", {
  style <- radial_style()

  expect_s3_class(style, "radial_style")
  expect_true(is.list(style))
  expect_equal(
    names(style),
    c(
      "bar_width", "bar_alpha", "bar_border_color", "category_label_size",
      "tick_label_size", "value_label_size", "spoke_color", "spoke_linetype",
      "grid_color", "strip_text_size", "radial_breaks", "radial_limits",
      "start_angle", "clockwise", "show_grid_circles", "show_value_labels",
      "legend_position"
    )
  )
  expect_equal(style$legend_position, "bottom")
})

test_that("radial_style() validates ranges and required fields", {
  expect_error(radial_style(bar_width = 0.09), "bar_width")
  expect_error(radial_style(bar_width = 1.01), "bar_width")
  expect_error(radial_style(bar_alpha = -0.01), "bar_alpha")
  expect_error(radial_style(bar_alpha = 1.01), "bar_alpha")
  expect_error(radial_style(show_grid_circles = 1), "show_grid_circles")
  expect_error(radial_style(radial_breaks = c(1, 2)), "include 0")
  expect_error(radial_style(radial_breaks = 0), "length >= 2")
  expect_error(radial_style(radial_limits = c(1, 1)), "less than")
  expect_error(radial_style(start_angle = NA_real_), "start_angle")
})

test_that("radial_plot() builds with default continuous fill and bare columns", {
  df <- data.frame(
    category = LETTERS[1:6],
    value = c(2, 4, 6, 8, 10, 12)
  )

  p <- radial_plot(df, category, value)

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$coordinates, "CoordPolar"))
  expect_equal(p$coordinates$clip, "off")
  expect_true(inherits(fill_scale(p), "ScaleContinuous"))
  expect_gt(y_scale(p)$limits[2], max(df$value))
  expect_true(all(c("GeomSegment", "GeomCol", "GeomLabel", "GeomText") %in% layer_geom_classes(p)))
  expect_gte(length(p$layers), 4L)
  expect_true(inherits(p$theme$axis.text, "element_blank"))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("radial_plot() accepts character column names for numeric fill", {
  df <- data.frame(
    category = LETTERS[1:4],
    value = c(2, 5, 7, 9),
    score = c(10, 30, 60, 90)
  )

  p <- radial_plot(
    df,
    category = "category",
    value = "value",
    fill = "score",
    palette = c("navy", "white", "firebrick"),
    limits = c(0, 100)
  )

  scale <- fill_scale(p)

  expect_s3_class(p, "ggplot")
  expect_true(inherits(scale, "ScaleContinuous"))
  expect_equal(scale$name, "score")
  expect_equal(scale$limits, c(0, 100))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("radial_plot() supports discrete fill palettes and validates mismatches", {
  df <- data.frame(
    category = LETTERS[1:4],
    value = c(3, 6, 9, 12),
    risk = factor(
      c("low", "medium", "high", "critical"),
      levels = c("low", "medium", "high", "critical")
    )
  )

  p_default <- radial_plot(df, category, value, fill = risk)
  p_named <- radial_plot(
    df,
    category,
    value,
    fill = risk,
    palette = c(
      low = "#2E7D32",
      medium = "#F0E442",
      high = "#E69F00",
      critical = "#C00000"
    )
  )
  p_unnamed <- radial_plot(
    df,
    category,
    value,
    fill = risk,
    palette = c("#2E7D32", "#F0E442", "#E69F00", "#C00000")
  )

  expect_true(inherits(fill_scale(p_default), "ScaleDiscrete"))
  expect_true(inherits(fill_scale(p_named), "ScaleDiscrete"))
  expect_true(inherits(fill_scale(p_unnamed), "ScaleDiscrete"))
  expect_equal(fill_scale(p_named)$name, "risk")
  expect_equal(fill_scale(p_named)$limits, levels(df$risk))
  expect_error(
    radial_plot(df, category, value, fill = risk, palette = "risk4"),
    "Unnamed `palette` must have length 4"
  )
  expect_error(
    radial_plot(
      df,
      category,
      value,
      fill = risk,
      palette = c(low = "#2E7D32", medium = "#F0E442", high = "#E69F00")
    ),
    "missing colors for level"
  )
  expect_error(
    radial_plot(df, category, value, fill = risk, palette = c("#2E7D32", "#F0E442")),
    "Unnamed `palette` must have length 4"
  )
})

test_that("radial_plot() facets cleanly and custom style settings are applied", {
  df <- data.frame(
    category = rep(LETTERS[1:6], 2),
    value = c(2, 4, 6, 8, 10, 12, 3, 5, 7, 9, 11, 13),
    scenario = rep(c("baseline", "stress"), each = 6)
  )
  style <- radial_style(
    show_value_labels = TRUE,
    show_grid_circles = FALSE,
    legend_position = "right"
  )

  p <- radial_plot(df, category, value, facet = scenario, style = style, facet_nrow = 1)

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
  expect_true(inherits(fill_scale(p), "ScaleContinuous"))
  expect_equal(p$theme$legend.position, "right")
  expect_true(inherits(p$theme$panel.grid.major.y, "element_blank"))
  expect_equal(sum(layer_geom_classes(p) == "GeomText"), 2L)
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("radial_plot() handles value NAs according to na_handling and drops NA-value rows", {
  df <- data.frame(
    category = c("A", NA, "C", "D"),
    value = c(4, NA_real_, 8, 10),
    scenario = c("S1", "S1", NA, "S2"),
    risk = c("low", NA, "high", "medium")
  )

  expect_silent(
    p_silent <- radial_plot(
      df,
      category,
      value,
      facet = scenario,
      fill = risk,
      na_handling = "silent"
    )
  )
  expect_warning(
    p_warn <- radial_plot(
      df,
      category,
      value,
      facet = scenario,
      fill = risk,
      na_handling = "warn"
    ),
    "contains 1 NA"
  )
  expect_error(
    radial_plot(df, category, value, facet = scenario, fill = risk, na_handling = "error"),
    "contains 1 NA"
  )

  expect_s3_class(p_silent, "ggplot")
  expect_s3_class(p_warn, "ggplot")
  expect_equal(nrow(p_silent$layers[[2]]$data), 3L)
  expect_false(anyNA(p_silent$layers[[2]]$data$.val))
  expect_true(inherits(fill_scale(p_silent), "ScaleDiscrete"))
})

test_that("radial_plot() validates empty data, non-numeric values, and unsupported fill types", {
  expect_error(
    radial_plot(data.frame(), category, value),
    "non-empty data.frame"
  )

  bad_value_df <- data.frame(category = c("A", "B"), value = c("1", "2"))
  expect_error(
    radial_plot(bad_value_df, category, value),
    "value` column must be numeric"
  )

  bad_fill_df <- data.frame(
    category = c("A", "B"),
    value = c(1, 2),
    fill_col = I(list(1, 2))
  )
  expect_error(
    radial_plot(bad_fill_df, category, value, fill = fill_col),
    "fill` column"
  )
})

test_that("radial_plot() edge cases around degenerate continuous ranges are explicit", {
  zero_df <- data.frame(category = LETTERS[1:3], value = c(0, 0, 0))

  expect_error(
    radial_plot(zero_df, category, value),
    "non-degenerate range"
  )

  p_zero <- radial_plot(zero_df, category, value, limits = c(0, 1))

  expect_s3_class(p_zero, "ggplot")
  expect_silent(ggplot2::ggplot_build(p_zero))

  one_row_df <- data.frame(category = "A", value = 5, band = "only")

  expect_error(
    radial_plot(one_row_df, category, value),
    "non-degenerate range"
  )

  p_one <- radial_plot(one_row_df, category, value, fill = band)

  expect_s3_class(p_one, "ggplot")
  expect_silent(ggplot2::ggplot_build(p_one))
})
