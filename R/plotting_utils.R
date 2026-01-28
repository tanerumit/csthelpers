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
#'   numeric column for stream order (see `stream_order_col`).
#' @param points_sf `sf` object (points). Optional point layer. Must contain a
#'   column for labels (see `label_col`).
#' @param stream_order_col Character scalar. Name of the numeric stream order
#'   column in `rivers_sf`. Default is `"strord"`.
#' @param label_col Character scalar. Name of the label column in `points_sf`.
#'   Default is `"name"`.
#' @param main_stem_order Numeric scalar. Minimum stream order threshold for
#'   filtering rivers. Rows with stream order >= this value are plotted
#'   (inclusive).
#' @param add_scale Logical. If `TRUE`, adds `ggspatial::annotation_scale()`.
#' @param add_north Logical. If `TRUE`, adds `ggspatial::annotation_north_arrow()`.
#' @param show_river_legend Logical. If `TRUE`, displays a legend for river
#'   line widths (stream order). Default is `FALSE`.
#'
#' @return A `ggplot` object with an additional attribute `"bbox"` containing
#'
#' @examples
#' \dontrun{
#' p <- plot_basin_map(basins_sf, rivers_sf, points_sf, main_stem_order = 2)
#' print(p)
#'
#' # Access the bounding box
#' attr(p, "bbox")
#' }
#'
#' @importFrom rlang .data
#' @importFrom sf st_crs st_transform st_coordinates st_bbox
#' @importFrom dplyr filter
#' @import ggplot2
#' @export
plot_basin_map <- function(
    basins_sf = NULL,
    rivers_sf = NULL,
    points_sf = NULL,
    stream_order_col = "strord",
    label_col = "name",
    main_stem_order = 2,
    add_scale = FALSE,
    add_north = FALSE,
    show_river_legend = FALSE
) {
  # ---------------------------------------------------------------------------

  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (is.null(basins_sf) && is.null(rivers_sf) && is.null(points_sf)) {
    stop(
      "At least one of `basins_sf`, `rivers_sf`, or `points_sf` must be provided.",
      call. = FALSE
    )
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
  if (!is.character(stream_order_col) || length(stream_order_col) != 1L || is.na(stream_order_col)) {
    stop("`stream_order_col` must be a non-missing character scalar.", call. = FALSE)
  }
  if (!is.character(label_col) || length(label_col) != 1L || is.na(label_col)) {
    stop("`label_col` must be a non-missing character scalar.", call. = FALSE)
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
  if (!is.logical(show_river_legend) || length(show_river_legend) != 1L || is.na(show_river_legend)) {
    stop("`show_river_legend` must be a non-missing logical scalar.", call. = FALSE)
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
  # CRS harmonization (use == for more forgiving comparison)
  # ---------------------------------------------------------------------------
  if (!is.null(basins_sf) && sf::st_crs(basins_sf) != ref_crs) {
    basins_sf <- sf::st_transform(basins_sf, ref_crs)
  }
  if (!is.null(rivers_sf) && sf::st_crs(rivers_sf) != ref_crs) {
    rivers_sf <- sf::st_transform(rivers_sf, ref_crs)
  }
  if (!is.null(points_sf) && sf::st_crs(points_sf) != ref_crs) {
    points_sf <- sf::st_transform(points_sf, ref_crs)
  }

  # ---------------------------------------------------------------------------
  # Compute combined bounding box for all layers
  # ---------------------------------------------------------------------------
  bbox_list <- list()
  if (!is.null(basins_sf)) bbox_list$basins <- sf::st_bbox(basins_sf)
  if (!is.null(rivers_sf)) bbox_list$rivers <- sf::st_bbox(rivers_sf)
  if (!is.null(points_sf)) bbox_list$points <- sf::st_bbox(points_sf)

  combined_bbox <- Reduce(
    f = function(a, b) {
      c(
        xmin = min(a[["xmin"]], b[["xmin"]]),
        ymin = min(a[["ymin"]], b[["ymin"]]),
        xmax = max(a[["xmax"]], b[["xmax"]]),
        ymax = max(a[["ymax"]], b[["ymax"]])
      )
    },
    x = bbox_list
  )
  class(combined_bbox) <- "bbox"
  attr(combined_bbox, "crs") <- ref_crs

  # ---------------------------------------------------------------------------
  # Prepare point coordinates + label nudges (based on overall map extent)
  # ---------------------------------------------------------------------------
  points_df <- NULL
  nudge_x <- 0
  nudge_y <- 0

  if (!is.null(points_sf)) {
    if (!(label_col %in% names(points_sf))) {
      stop(
        sprintf("`points_sf` must contain a `%s` column for labelling.", label_col),
        call. = FALSE
      )
    }

    coords <- sf::st_coordinates(points_sf)
    points_df <- cbind(sf::st_drop_geometry(points_sf), coords)
    # Rename label column to internal name for consistent aes mapping
    points_df[[".label"]] <- points_df[[label_col]]

    # Use the combined bounding box for nudge calculation (more stable).
    bbox_width  <- as.numeric(combined_bbox[["xmax"]] - combined_bbox[["xmin"]])
    bbox_height <- as.numeric(combined_bbox[["ymax"]] - combined_bbox[["ymin"]])

    # If degenerate extent (single point), fall back to small fixed nudges.
    NUDGE_X_FRAC <- 0.03
    NUDGE_Y_FRAC <- 0.04
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
  LABEL_FILL        <- ggplot2::alpha("white", 0.6)

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
    if (!(stream_order_col %in% names(rivers_sf))) {
      stop(
        sprintf(
          "`rivers_sf` must contain a `%s` column for stream order.",
          stream_order_col
        ),
        call. = FALSE
      )
    }
    if (!is.numeric(rivers_sf[[stream_order_col]])) {
      stop(
        sprintf("`rivers_sf$%s` must be numeric.", stream_order_col),
        call. = FALSE
      )
    }

    # Rename to internal column for consistent aes mapping
    rivers_sf[[".stream_order"]] <- rivers_sf[[stream_order_col]]
    rivers_filt <- dplyr::filter(rivers_sf, .data$.stream_order >= main_stem_order)

    if (nrow(rivers_filt) == 0) {
      warning(
        sprintf(
          "No rivers remain after filtering with `main_stem_order = %s`. Consider lowering the threshold.",
          main_stem_order
        ),
        call. = FALSE
      )
    } else {
      legend_guide <- if (show_river_legend) "legend" else "none"

      p <- p +
        geom_sf(
          data = rivers_filt,
          aes(linewidth = .data$.stream_order),
          color = RIVER_COLOR,
          alpha = RIVER_ALPHA
        ) +
        scale_linewidth(
          range = RIVER_LWD_RANGE,
          guide = legend_guide,
          name = "Stream order"
        )
    }
  }

  if (!is.null(points_df)) {
    p <- p +
      geom_label(
        data = points_df,
        aes(x = .data$X, y = .data$Y, label = .data$.label),
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
    if (!requireNamespace("ggspatial", quietly = TRUE)) {
      warning(
        "Package 'ggspatial' is required for `add_scale = TRUE` but is not installed.",
        call. = FALSE
      )
    } else {
      p <- p + ggspatial::annotation_scale(
        location = "bl",
        width_hint = 0.3,
        text_cex = 0.75
      )
    }
  }

  if (isTRUE(add_north)) {
    if (!requireNamespace("ggspatial", quietly = TRUE)) {
      warning(
        "Package 'ggspatial' is required for `add_north = TRUE` but is not installed.",
        call. = FALSE
      )
    } else {
      p <- p + ggspatial::annotation_north_arrow(
        location = "tl",
        which_north = "true",
        style = ggspatial::north_arrow_fancy_orienteering,
        height = grid::unit(1, "cm"),
        width  = grid::unit(1, "cm")
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Final layout (use expand = TRUE to avoid clipping labels)
  # ---------------------------------------------------------------------------
  p <- p +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    coord_sf(expand = TRUE)

  # Attach bounding box as attribute
  attr(p, "bbox") <- combined_bbox

  p
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
#' @param value_col Numeric variable specifying the radial magnitude (bar length).
#'   `NA` values are permitted and will result in no bar being drawn for that
#'   observation; a warning is issued if any `NA`s are present.
#' @param color_col Variable defining the categorical fill grouping.
#' @param cat_colors Named character vector defining fill colors for `color_col` levels.
#'   Names must match the factor levels to enforce consistent ordering.
#' @param radial_breaks Optional numeric vector of radial breaks (must include `0`).
#'   If `NULL`, breaks are computed via `pretty()`.
#' @param start_angle Numeric. Angular position (degrees) of the `"GAP"` sector where
#'   tick labels are printed. Default is `90` (top). Values are normalized to `[0, 360)`.
#' @param direction Numeric. `1` for counter-clockwise (default), `-1` for clockwise.
#' @param facet_nrow Integer. Number of rows for facet layout. Default `NULL` allows
#'   `facet_wrap()` to choose automatically.
#' @param legend_title Character. Title for the fill legend. Default is `"Category"`.
#' @param spoke_to_max Logical. If `TRUE` (default), dashed helper spokes extend to
#'   `max(radial_breaks)`. If `FALSE`, spokes extend only to each bar's value.
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
#' cat_colors <- c("No" = "green", "Low" = "yellow", "Medium" = "orange", "High" = "red")
#'
#' sos_radialplot(
#'   data = df,
#'   axis_col = group,
#'   panel_col = scenario,
#'   value_col = value,
#'   color_col = risk_class,
#'   cat_colors = cat_colors,
#'   radial_breaks = seq(0, 50, 10),
#'   start_angle = 90,
#'   legend_title = "Risk Class"
#' )
#' }
#'
#' @importFrom rlang .data
#' @importFrom grid unit
#' @importFrom dplyr rename group_by slice ungroup mutate bind_rows filter
#' @import ggplot2
#' @export
sos_radialplot <- function(
    data,
    axis_col,
    panel_col,
    value_col,
    color_col,
    cat_colors,
    radial_breaks = NULL,
    start_angle = 90,
    direction = 1,
    facet_nrow = NULL,
    legend_title = "Category",
    spoke_to_max = TRUE
) {

  # ---------------------------------------------------------------------------
  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` must have at least one row.", call. = FALSE)
  }
  if (!is.numeric(start_angle) || length(start_angle) != 1L || is.na(start_angle)) {
    stop("`start_angle` must be a non-missing numeric scalar (degrees).", call. = FALSE)
  }
  if (!is.numeric(direction) || length(direction) != 1L || !direction %in% c(-1, 1)) {
    stop("`direction` must be either 1 (counter-clockwise) or -1 (clockwise).", call. = FALSE)
  }
  if (!is.null(facet_nrow) && (!is.numeric(facet_nrow) || length(facet_nrow) != 1L ||
                               is.na(facet_nrow) || facet_nrow < 1L || facet_nrow != as.integer(facet_nrow))) {
    stop("`facet_nrow` must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.character(legend_title) || length(legend_title) != 1L || is.na(legend_title)) {
    stop("`legend_title` must be a non-missing character scalar.", call. = FALSE)
  }
  if (!is.logical(spoke_to_max) || length(spoke_to_max) != 1L || is.na(spoke_to_max)) {
    stop("`spoke_to_max` must be a non-missing logical scalar.", call. = FALSE)
  }
  if (!is.character(cat_colors) || is.null(names(cat_colors)) || anyNA(names(cat_colors))) {
    stop("`cat_colors` must be a *named* character vector (names = category levels).", call. = FALSE)
  }

  # Normalize start_angle to [0, 360)
  start_angle <- start_angle %% 360

  # ---------------------------------------------------------------------------
  # Standardize column names (tidy-eval) for internal plotting
  # Use distinctive internal names to avoid collision with user data
  # ---------------------------------------------------------------------------
  df <- dplyr::rename(
    data,
    .axis  = {{ axis_col }},
    .panel = {{ panel_col }},
    .value = {{ value_col }},
    .color = {{ color_col }}
  )

  required_cols <- c(".axis", ".panel", ".value", ".color")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      "Column mapping failed. Could not find columns for: ",
      paste(gsub("^\\.", "", missing_cols), collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(df[[".value"]])) {
    stop("`value_col` must be numeric.", call. = FALSE)
  }

  # Check for all-NA values
  if (all(is.na(df[[".value"]]))) {
    stop("`value_col` contains only NA values.", call. = FALSE)
  }

  # Warn about NA values
  n_na <- sum(is.na(df[[".value"]]))
  if (n_na > 0L) {
    warning(
      sprintf("`value_col` contains %d NA value(s); these observations will not be plotted.", n_na),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Factor ordering: preserve input order for axis/panels; enforce cat_colors for fill
  # Use ordered = FALSE to avoid unexpected scale/stat interactions
  # ---------------------------------------------------------------------------
  df[[".axis"]]  <- factor(df[[".axis"]],  levels = unique(df[[".axis"]]))
  df[[".panel"]] <- factor(df[[".panel"]], levels = unique(df[[".panel"]]))

  # Ensure all observed color levels exist in cat_colors; enforce consistent ordering.
  observed_levels <- unique(as.character(df[[".color"]]))
  observed_levels <- observed_levels[!is.na(observed_levels)]
  unknown_levels <- setdiff(observed_levels, names(cat_colors))
  if (length(unknown_levels) > 0L) {
    stop(
      "`cat_colors` is missing names for these `color_col` levels: ",
      paste(unknown_levels, collapse = ", "),
      call. = FALSE
    )
  }
  df[[".color"]] <- factor(df[[".color"]], levels = names(cat_colors))

  axis_levels  <- levels(df[[".axis"]])
  panel_levels <- levels(df[[".panel"]])

  # ---------------------------------------------------------------------------
  # Radial breaks: compute before any tick data are built
  # ---------------------------------------------------------------------------
  if (is.null(radial_breaks)) {
    max_value <- max(df[[".value"]], na.rm = TRUE)
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
  # Explicitly set .color = NA to avoid inheriting misleading values
  # ---------------------------------------------------------------------------
  gap_rows <- df |>
    dplyr::group_by(.data$.panel) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::mutate(.axis = "GAP", .value = 0, .color = NA)

  df <- dplyr::bind_rows(df, gap_rows)
  df[[".axis"]] <- factor(df[[".axis"]], levels = c(axis_levels, "GAP"))

  # ---------------------------------------------------------------------------
  # Label data for axes around rim (exclude GAP as a label)
  # ---------------------------------------------------------------------------
  lab_df <- data.frame(
    .axis = factor(axis_levels, levels = c(axis_levels, "GAP"))
  )

  # Tick labels placed at the GAP sector for each panel.
  tick_df <- base::expand.grid(
    .panel = panel_levels,
    .axis  = factor("GAP", levels = c(axis_levels, "GAP")),
    .y     = radial_breaks[radial_breaks != 0],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Prepare spoke data
  # ---------------------------------------------------------------------------
  spoke_df <- dplyr::filter(df, .data$.axis != "GAP", !is.na(.data$.value))

  if (spoke_to_max) {
    spoke_df[[".yend"]] <- ymax_ext
  } else {
    spoke_df[[".yend"]] <- spoke_df[[".value"]]
  }

  # ---------------------------------------------------------------------------
  # Theme + styling constants
  # ---------------------------------------------------------------------------
  LABEL_SIZE <- 3

  ggtheme <- theme_bw() +
    theme(
      axis.text.x       = element_blank(),
      axis.ticks        = element_blank(),
      axis.text.y       = element_blank(),
      axis.title        = element_text(size = 12, color = "black"),
      panel.spacing     = grid::unit(1, "lines"),
      plot.title.position = "plot",
      # Add margin to prevent label clipping with clip = "off"
      plot.margin       = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")
    )

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot(df) +
    ggtheme +

    # Dashed helper spokes (exclude GAP)
    geom_segment(
      data = spoke_df,
      aes(x = .data$.axis, xend = .data$.axis, y = 0, yend = .data$.yend),
      linetype = "dashed",
      color = "gray12"
    ) +

    # Bars (exclude GAP and NA values)
    geom_col(
      data = dplyr::filter(df, .data$.axis != "GAP", !is.na(.data$.value)),
      aes(x = .data$.axis, y = .data$.value, fill = .data$.color),
      width = 0.8,
      color = "gray20",
      alpha = 0.9,
      show.legend = TRUE
    ) +

    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_continuous(limits = c(0, ymax_ext), breaks = radial_breaks) +
    scale_fill_manual(values = cat_colors, name = legend_title, drop = FALSE, na.translate = FALSE) +

    facet_wrap(~ .panel, nrow = facet_nrow) +
    coord_polar(clip = "off", start = start_angle * pi / 180, direction = direction) +
    labs(x = NULL, y = NULL) +

    # Axis labels around rim (fixed radius)
    geom_label(
      data = lab_df,
      aes(x = .data$.axis, y = ymax_ext, label = .data$.axis),
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
      aes(x = .data$.axis, y = .data$.y, label = .data$.y),
      size = LABEL_SIZE,
      color = "black",
      alpha = 1,
      inherit.aes = FALSE
    )

  p
}
