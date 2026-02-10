test_that("plot_basin_map() fails if all inputs are NULL", {
  expect_error(
    plot_basin_map(),
    "At least one of `basins_sf`, `rivers_sf`, or `points_sf` must be provided."
  )
})

test_that("plot_basin_map() validates sf inputs types", {
  expect_error(plot_basin_map(basins_sf = data.frame(a = 1)), "`basins_sf` must be an sf object")
  expect_error(plot_basin_map(rivers_sf = data.frame(a = 1)), "`rivers_sf` must be an sf object")
  expect_error(plot_basin_map(points_sf = data.frame(a = 1)), "`points_sf` must be an sf object")
})

test_that("plot_basin_map() validates scalar argument types", {
  basins <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
    ))), crs = 4326)
  )

  expect_error(plot_basin_map(basins_sf = basins, stream_order_col = NA_character_), "stream_order_col")
  expect_error(plot_basin_map(basins_sf = basins, stream_order_col = c("a", "b")), "stream_order_col")
  expect_error(plot_basin_map(basins_sf = basins, label_col = NA_character_), "label_col")
  expect_error(plot_basin_map(basins_sf = basins, label_col = c("a", "b")), "label_col")

  expect_error(plot_basin_map(basins_sf = basins, main_stem_order = NA_real_), "main_stem_order")
  expect_error(plot_basin_map(basins_sf = basins, main_stem_order = c(1, 2)), "main_stem_order")
  expect_error(plot_basin_map(basins_sf = basins, main_stem_order = "2"), "main_stem_order")

  expect_error(plot_basin_map(basins_sf = basins, add_scale = NA), "add_scale")
  expect_error(plot_basin_map(basins_sf = basins, add_scale = 1), "add_scale")

  expect_error(plot_basin_map(basins_sf = basins, add_north = NA), "add_north")
  expect_error(plot_basin_map(basins_sf = basins, add_north = 1), "add_north")

  expect_error(plot_basin_map(basins_sf = basins, show_river_legend = NA), "show_river_legend")
  expect_error(plot_basin_map(basins_sf = basins, show_river_legend = 1), "show_river_legend")
})

test_that("plot_basin_map() errors if CRS is missing/NA on the first non-NULL layer", {
  basins_geom <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
  ))), crs = 4326)
  sf::st_crs(basins_geom) <- NA
  basins_na_crs <- sf::st_sf(id = 1, geometry = basins_geom)

  expect_error(
    plot_basin_map(basins_sf = basins_na_crs),
    "Could not determine a valid CRS from inputs \\(missing/NA CRS\\)"
  )
})

test_that("plot_basin_map() returns ggplot with bbox attribute", {
  basins <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(2, 0), c(2, 1), c(0, 1), c(0, 0)
    ))), crs = 4326)
  )

  p <- plot_basin_map(basins_sf = basins)
  expect_s3_class(p, "ggplot")

  bb <- attr(p, "bbox")
  expect_true(inherits(bb, "bbox"))
  expect_true(all(c("xmin", "ymin", "xmax", "ymax") %in% names(bb)))
  expect_true(!is.null(attr(bb, "crs")))
  expect_true(sf::st_crs(attr(bb, "crs")) == sf::st_crs(basins))
})

test_that("plot_basin_map() computes combined bbox across layers", {
  basins <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
    ))), crs = 4326)
  )

  rivers <- sf::st_sf(
    strord = c(1, 3),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(-10, 0), c(-9, 1))),
      sf::st_linestring(rbind(c(5, 0), c(6, 1))),
      crs = 4326
    )
  )

  points <- sf::st_sf(
    name = "P1",
    geometry = sf::st_sfc(sf::st_point(c(100, 100)), crs = 4326)
  )

  p <- plot_basin_map(basins_sf = basins, rivers_sf = rivers, points_sf = points)
  bb <- attr(p, "bbox")

  expect_equal(bb[["xmin"]], -10)
  expect_equal(bb[["ymin"]], 0)
  expect_equal(bb[["xmax"]], 100)
  expect_equal(bb[["ymax"]], 100)
})

test_that("plot_basin_map() transforms CRS to reference CRS (first non-NULL layer)", {
  basins_4326 <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
    ))), crs = 4326)
  )

  # points in Web Mercator (meters)
  points_3857 <- sf::st_sf(
    name = "P1",
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 3857)
  )

  p <- plot_basin_map(basins_sf = basins_4326, points_sf = points_3857)
  bb <- attr(p, "bbox")
  expect_true(sf::st_crs(attr(bb, "crs")) == sf::st_crs(basins_4326))
})

test_that("plot_basin_map() errors when points layer is missing label column", {
  points <- sf::st_sf(
    wrong = "P1",
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326)
  )

  expect_error(
    plot_basin_map(points_sf = points, label_col = "name"),
    "`points_sf` must contain a `name` column for labelling"
  )
})

test_that("plot_basin_map() errors when rivers layer is missing stream order column", {
  rivers <- sf::st_sf(
    wrong = 1,
    geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(1, 1))), crs = 4326)
  )

  expect_error(
    plot_basin_map(rivers_sf = rivers, stream_order_col = "strord"),
    "`rivers_sf` must contain a `strord` column for stream order"
  )
})

test_that("plot_basin_map() errors when stream order column is non-numeric", {
  rivers <- sf::st_sf(
    strord = c("1", "2"),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(0, 0), c(1, 1))),
      sf::st_linestring(rbind(c(1, 0), c(2, 1))),
      crs = 4326
    )
  )

  expect_error(
    plot_basin_map(rivers_sf = rivers),
    "`rivers_sf\\$strord` must be numeric"
  )
})

test_that("plot_basin_map() warns when no rivers remain after filtering", {
  rivers <- sf::st_sf(
    strord = c(1, 1),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(0, 0), c(1, 1))),
      sf::st_linestring(rbind(c(1, 0), c(2, 1))),
      crs = 4326
    )
  )

  expect_warning(
    plot_basin_map(rivers_sf = rivers, main_stem_order = 10),
    "No rivers remain after filtering"
  )
})

test_that("plot_basin_map() does not error for degenerate point extent (single point)", {
  points <- sf::st_sf(
    name = "P1",
    geometry = sf::st_sfc(sf::st_point(c(5, 5)), crs = 4326)
  )

  p <- plot_basin_map(points_sf = points)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(attr(p, "bbox"), "bbox"))
})

test_that("plot_basin_map() adds scale/north warnings if ggspatial is missing", {
  # Only run these expectations if ggspatial is NOT installed, otherwise skip.
  if (requireNamespace("ggspatial", quietly = TRUE)) {
    skip("ggspatial installed; missing-package warnings are not expected.")
  }

  basins <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
    ))), crs = 4326)
  )

  expect_warning(
    plot_basin_map(basins_sf = basins, add_scale = TRUE),
    "Package 'ggspatial' is required for `add_scale = TRUE`"
  )

  expect_warning(
    plot_basin_map(basins_sf = basins, add_north = TRUE),
    "Package 'ggspatial' is required for `add_north = TRUE`"
  )
})
