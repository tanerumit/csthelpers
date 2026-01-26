#' Plot a basin overview map (basins, river network, and labelled points)
#'
#' @description
#' Creates a simple basin overview map using `geom_sf()`, with optional
#' basin polygons, river network lines filtered by stream order, and labelled
#' point locations.
#'
#' The function harmonizes CRS across inputs by transforming all layers to the
#' CRS of the first non-`NULL` object among `basins_sf`, `rivers_sf`, `points_sf`.
#'
#' @param basins_sf `sf` object (polygons). Optional basin boundaries.
#' @param rivers_sf `sf` object (lines). Optional river network. Must contain a
#'   numeric column named `strord` (stream order).
#' @param points_sf `sf` object (points). Optional point layer. Must contain a
#'   column named `name` for labels.
#' @param main_stem_order Numeric scalar. Minimum stream order threshold used
#'   for filtering rivers: rows with `strord > main_stem_order` are plotted.
#' @param add_scale Logical. If `TRUE`, adds `ggspatial::annotation_scale()`.
#' @param add_north Logical. If `TRUE`, adds `ggspatial::annotation_north_arrow()`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' p <- plot_basin_map(basins_sf, rivers_sf, points_sf, main_stem_order = 2)
#' print(p)
#' }
#'
#' @import ggplot2 dplyr
#' @export
plot_basin_map <- function(
    basins_sf = NULL,
    rivers_sf = NULL,
    points_sf = NULL,
    main_stem_order = 2,
    add_scale = FALSE,
    add_north = FALSE
) {
  # ---------------------------------------------------------------------------
  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (is.null(basins_sf) && is.null(rivers_sf) && is.null(points_sf)) {
    stop("At least one of `basins_sf`, `rivers_sf`, or `points_sf` must be provided.", call. = FALSE)
  }
  if (!is.null(basins_sf) && !inherits(basins_sf, "sf")) {
    stop("`basins_sf` must be an sf object.", call. = FALSE)
  }
  if (!is.null(rivers_sf) && !inherits(rivers_sf, "sf")) {
    stop("`rivers_sf` must be an sf object.", call. = FALSE)
  }
  if (!is.null(points_sf) && !inherits(points_sf, "sf")) {
    stop("`points_sf` must be an sf object.", call. = FALSE)
  }
  if (!is.numeric(main_stem_order) || length(main_stem_order) != 1L || is.na(main_stem_order)) {
    stop("`main_stem_order` must be a non-missing numeric scalar.", call. = FALSE)
  }
  if (!is.logical(add_scale) || length(add_scale) != 1L || is.na(add_scale)) {
    stop("`add_scale` must be a non-missing logical scalar.", call. = FALSE)
  }
  if (!is.logical(add_north) || length(add_north) != 1L || is.na(add_north)) {
    stop("`add_north` must be a non-missing logical scalar.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Determine reference CRS
  # ---------------------------------------------------------------------------
  ref_crs <- NULL
  if (!is.null(basins_sf)) {
    ref_crs <- sf::st_crs(basins_sf)
  } else if (!is.null(rivers_sf)) {
    ref_crs <- sf::st_crs(rivers_sf)
  } else if (!is.null(points_sf)) {
    ref_crs <- sf::st_crs(points_sf)
  }

  if (is.null(ref_crs) || is.na(ref_crs)) {
    stop("Could not determine a valid CRS from inputs (missing/NA CRS).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # CRS harmonization
  # ---------------------------------------------------------------------------
  if (!is.null(basins_sf) && !identical(sf::st_crs(basins_sf), ref_crs)) {
    basins_sf <- sf::st_transform(basins_sf, ref_crs)
  }
  if (!is.null(rivers_sf) && !identical(sf::st_crs(rivers_sf), ref_crs)) {
    rivers_sf <- sf::st_transform(rivers_sf, ref_crs)
  }
  if (!is.null(points_sf) && !identical(sf::st_crs(points_sf), ref_crs)) {
    points_sf <- sf::st_transform(points_sf, ref_crs)
  }

  # ---------------------------------------------------------------------------
  # Prepare point coordinates + label nudges (based on overall map extent)
  # ---------------------------------------------------------------------------
  points_df <- NULL
  nudge_x <- 0
  nudge_y <- 0

  if (!is.null(points_sf)) {
    if (!("name" %in% names(points_sf))) {
      stop("`points_sf` must contain a `name` column for labelling.", call. = FALSE)
    }

    coords <- sf::st_coordinates(points_sf)
    points_df <- cbind(points_sf, coords)

    # Use the map bounding box (more stable than point-to-point range).
    bbox <- sf::st_bbox(points_sf)
    bbox_width  <- as.numeric(bbox[["xmax"]] - bbox[["xmin"]])
    bbox_height <- as.numeric(bbox[["ymax"]] - bbox[["ymin"]])

    # If only one point (or degenerate extent), fall back to small fixed nudges.
    NUDGE_X_FRAC <- 0.05
    NUDGE_Y_FRAC <- 0.07
    nudge_x <- if (is.finite(bbox_width) && bbox_width > 0) bbox_width * NUDGE_X_FRAC else 0
    nudge_y <- if (is.finite(bbox_height) && bbox_height > 0) bbox_height * NUDGE_Y_FRAC else 0
  }

  # ---------------------------------------------------------------------------
  # Plot styling constants
  # ---------------------------------------------------------------------------
  BASIN_FILL        <- "grey95"
  BASIN_BORDER      <- "grey70"
  BASIN_LINEWIDTH   <- 0.2

  RIVER_COLOR       <- "steelblue4"
  RIVER_ALPHA       <- 0.9
  RIVER_LWD_RANGE   <- c(0.2, 1.5)

  LABEL_TEXT_SIZE   <- 3
  LABEL_FILL        <- scales::alpha("white", 0.6)

  POINT_SHAPE       <- 21
  POINT_FILL        <- "black"
  POINT_BORDER      <- "white"
  POINT_SIZE        <- 2.5
  POINT_STROKE      <- 0.3

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot()

  if (!is.null(basins_sf)) {
    p <- p + geom_sf(
      data = basins_sf,
      fill = BASIN_FILL,
      color = BASIN_BORDER,
      linewidth = BASIN_LINEWIDTH
    )
  }

  if (!is.null(rivers_sf)) {
    if (!("strord" %in% names(rivers_sf))) {
      stop("`rivers_sf` must contain a numeric `strord` column for stream order.", call. = FALSE)
    }
    if (!is.numeric(rivers_sf[["strord"]])) {
      stop("`rivers_sf$strord` must be numeric.", call. = FALSE)
    }

    rivers_filt <- filter(rivers_sf, .data$strord > main_stem_order)

    p <- p +
      geom_sf(
        data = rivers_filt,
        aes(linewidth = .data$strord),
        color = RIVER_COLOR,
        alpha = RIVER_ALPHA
      ) +
      scale_linewidth(range = RIVER_LWD_RANGE, guide = "none")
  }

  if (!is.null(points_df)) {
    p <- p +
      geom_label(
        data = points_df,
        aes(x = .data$X, y = .data$Y, label = .data$name),
        size = LABEL_TEXT_SIZE,
        label.size = 0,
        fill = LABEL_FILL,
        color = "black",
        label.r = grid::unit(0.1, "lines"),
        nudge_y = nudge_y,
        nudge_x = nudge_x
      ) +
      geom_point(
        data = points_df,
        aes(x = .data$X, y = .data$Y),
        shape = POINT_SHAPE,
        fill = POINT_FILL,
        color = POINT_BORDER,
        size = POINT_SIZE,
        stroke = POINT_STROKE
      )
  }

  # ---------------------------------------------------------------------------
  # Optional map annotations
  # ---------------------------------------------------------------------------
  if (isTRUE(add_scale)) {
    p <- p + ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.3,
      text_cex = 0.75
    )
  }

  if (isTRUE(add_north)) {
    p <- p + ggspatial::annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering,
      height = grid::unit(1, "cm"),
      width  = grid::unit(1, "cm")
    )
  }

  # ---------------------------------------------------------------------------
  # Final layout
  # ---------------------------------------------------------------------------
  p +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    coord_sf(expand = FALSE)
}


#' Create a facetted radial (spider) plot with fixed tick placement
#'
#' @description
#' `sos_radialplot()` generates a facetted radial (polar) bar plot for visualizing
#' categorical groups arranged around a circle (e.g., indicators, sites, or components)
#' and their corresponding numeric values under different panels (e.g., scenarios).
#'
#' The function enforces consistent radial tick placement by reserving a dummy sector
#' (`"GAP"`) at a fixed angular position. The `"GAP"` sector is not drawn as a bar/spoke,
#' and is used only for printing tick labels.
#'
#' @param data A data frame containing the variables to plot.
#' @param axis_col Variable mapped to the circular axis (e.g., `group` or `site`).
#' @param panel_col Variable used to facet panels (e.g., `scenario`).
#' @param panel_var Optional variable for nested paneling. Currently unused; retained
#'   for interface stability.
#' @param value_col Numeric variable specifying the radial magnitude (bar length).
#' @param cat_colors Named character vector defining fill colors for `color_col` levels.
#'   Names must match the factor levels to enforce consistent ordering.
#' @param color_col Variable defining the categorical fill grouping.
#' @param radial_breaks Optional numeric vector of radial breaks (must include `0`).
#'   If `NULL`, breaks are computed via `pretty()`.
#' @param start_angle Numeric. Angular position (degrees) of the `"GAP"` sector where
#'   tick labels are printed. Default is `90` (top).
#'
#' @details
#' Internally, a dummy `"GAP"` level is appended to the axis factor levels to reserve
#' a stable angular position for tick labels. This makes facet-to-facet comparison
#' reliable even when the number of axis categories varies.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   group = rep(c("HF", "MF", "LF", "RFC", "IF", "Overall"), times = 4),
#'   scenario = rep(c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"), each = 6),
#'   value = c(5, 10, 25, 0, 5, 10,
#'             8, 12, 30, 10, 0, 15,
#'             10, 18, 35, 15, 5, 20,
#'             15, 20, 40, 18, 10, 25),
#'   risk_class = rep(c("Low", "Low", "Medium", "No", "Low", "Low"), times = 4)
#' )
#'
#' cat_colors <- c("No"="green", "Low"="yellow", "Medium"="orange", "High"="red")
#'
#' sos_radialplot(
#'   data = df,
#'   axis_col = group,
#'   panel_col = scenario,
#'   value_col = value,
#'   color_col = risk_class,
#'   cat_colors = cat_colors,
#'   radial_breaks = seq(0, 50, 10),
#'   start_angle = 90
#' )
#' }
#'
#' @import ggplot2 dplyr tibble
#' @export
sos_radialplot <- function(
    data,
    axis_col,
    panel_col,
    panel_var = NULL,
    value_col,
    cat_colors,
    color_col,
    radial_breaks = NULL,
    start_angle = 90
) {
  # ---------------------------------------------------------------------------
  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.numeric(start_angle) || length(start_angle) != 1L || is.na(start_angle)) {
    stop("`start_angle` must be a non-missing numeric scalar (degrees).", call. = FALSE)
  }
  if (!is.character(cat_colors) || is.null(names(cat_colors)) || anyNA(names(cat_colors))) {
    stop("`cat_colors` must be a *named* character vector (names = category levels).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Standardize column names (tidy-eval) for internal plotting
  # ---------------------------------------------------------------------------
  df <- rename(
    data,
    axis_col  = {{ axis_col }},
    panel_col = {{ panel_col }},
    value_col = {{ value_col }},
    color_col = {{ color_col }}
  )

  required_cols <- c("axis_col", "panel_col", "value_col", "color_col")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns after rename: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(df[["value_col"]])) {
    stop("`value_col` must be numeric.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Factor ordering: preserve input order for axis/panels; enforce cat_colors for fill
  # ---------------------------------------------------------------------------
  df[["axis_col"]]  <- factor(df[["axis_col"]],  levels = unique(df[["axis_col"]]),  ordered = TRUE)
  df[["panel_col"]] <- factor(df[["panel_col"]], levels = unique(df[["panel_col"]]), ordered = TRUE)

  # Ensure all observed color levels exist in cat_colors; enforce consistent ordering.
  observed_levels <- unique(as.character(df[["color_col"]]))
  unknown_levels <- setdiff(observed_levels, names(cat_colors))
  if (length(unknown_levels) > 0L) {
    stop(
      "`cat_colors` is missing names for these `color_col` levels: ",
      paste(unknown_levels, collapse = ", "),
      call. = FALSE
    )
  }
  df[["color_col"]] <- factor(df[["color_col"]], levels = names(cat_colors), ordered = TRUE)

  axis_levels  <- levels(df[["axis_col"]])
  panel_levels <- levels(df[["panel_col"]])

  # ---------------------------------------------------------------------------
  # Radial breaks: compute before any tick data are built
  # ---------------------------------------------------------------------------
  if (is.null(radial_breaks)) {
    max_value <- max(df[["value_col"]], na.rm = TRUE)
    radial_breaks <- pretty(c(0, max_value))
  } else {
    if (!is.numeric(radial_breaks) || length(radial_breaks) < 2L || anyNA(radial_breaks)) {
      stop("`radial_breaks` must be a numeric vector of length >= 2 with no NA.", call. = FALSE)
    }
    if (!any(radial_breaks == 0)) {
      stop("`radial_breaks` must include 0.", call. = FALSE)
    }
    radial_breaks <- sort(unique(radial_breaks))
  }
  ymax_ext <- max(radial_breaks)

  # ---------------------------------------------------------------------------
  # Add a GAP row per panel (for tick placement only)
  # ---------------------------------------------------------------------------
  gap_row <- df |>
    group_by(.data$panel_col) |>
    slice(1) |>
    ungroup() |>
    mutate(axis_col = "GAP", value_col = 0)

  df <- bind_rows(df, gap_row)
  df[["axis_col"]] <- factor(df[["axis_col"]], levels = c(axis_levels, "GAP"), ordered = TRUE)

  # ---------------------------------------------------------------------------
  # Label data for axes around rim (exclude GAP as a label)
  # ---------------------------------------------------------------------------
  lab_df <- tibble::tibble(
    axis_col = factor(axis_levels, levels = c(axis_levels, "GAP"), ordered = TRUE)
  )

  # Tick labels placed at the GAP sector for each panel.
  # Use base::expand.grid to avoid adding tidyr as a dependency.
  tick_df <- base::expand.grid(
    panel_col = panel_levels,
    axis_col  = factor("GAP", levels = c(axis_levels, "GAP"), ordered = TRUE),
    y         = radial_breaks[radial_breaks != 0],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Theme + styling constants
  # ---------------------------------------------------------------------------
  LABEL_SIZE <- 3

  ggtheme <- theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.y = element_blank(),
      axis.title  = element_text(size = 12, color = "black"),
      panel.spacing = grid::unit(1, "lines"),
      plot.title.position = "plot"
    )

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot(df) +
    ggtheme +
    # Dashed helper spokes (exclude GAP)
    geom_segment(
      data = filter(df, .data$axis_col != "GAP"),
      aes(x = .data$axis_col, xend = .data$axis_col, y = 0, yend = ymax_ext),
      linetype = "dashed",
      color = "gray12"
    ) +

    # Bars (exclude GAP)
    geom_col(
      data = filter(df, .data$axis_col != "GAP"),
      aes(x = .data$axis_col, y = .data$value_col, fill = .data$color_col),
      width = 0.8,
      color = "gray20",
      alpha = 0.9,
      show.legend = TRUE
    ) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_continuous(limits = c(0, ymax_ext), breaks = radial_breaks) +
    scale_fill_manual(values = cat_colors, name = "Levels", drop = FALSE) +
    facet_wrap(~ panel_col, nrow = 2) +
    coord_polar(clip = "on", start = start_angle * pi / 180) +
    labs(x = NULL, y = NULL) +

    # Axis labels around rim (fixed radius)
    geom_label(
      data = lab_df,
      aes(x = .data$axis_col, y = ymax_ext, label = .data$axis_col),
      fill = "white",
      label.size = 0,
      color = "black",
      alpha = 0.8,
      fontface = "bold",
      size = LABEL_SIZE,
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    # Radial tick labels at GAP
    geom_text(
      data = tick_df,
      aes(x = .data$axis_col, y = .data$y, label = .data$y),
      size = LABEL_SIZE,
      color = "black",
      alpha = 1,
      inherit.aes = FALSE
    )

  p
}
