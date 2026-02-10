#' Plot a basin overview map (basins, river network, and labelled points)
#'
#' @description
#' Creates a simple basin overview map using `geom_sf()`, with optional
#' basin polygons, river network lines filtered by stream order, and labelled
#' point locations.
#'
#' The function harmonizes CRS across inputs by transforming all layers to the
#' CRS of the first non-`NULL` object among `basins_sf`, `subbasins_sf`,
#' `rivers_sf`, `points_sf`.
#'
#' @param basins_sf `sf` object (polygons). Optional basin boundaries.
#' @param subbasins_sf `sf` object (polygons). Optional subbasin boundaries to
#'   overlay on the map.
#' @param rivers_sf `sf` object (lines). Optional river network. Must contain a
#'   numeric column for stream order (see `stream_order_col`).
#' @param points_sf `sf` object (points). Optional point layer. Must contain a
#'   column for labels (see `label_col`).
#' @param stream_order_col Character scalar. Name of the numeric stream order
#'   column in `rivers_sf`. Default is `"strord"`.
#' @param label_col Character scalar. Name of the label column in `points_sf`.
#'   Default is `"name"`.
#' @param min_stream_order Numeric scalar. Minimum stream order threshold for
#'   filtering rivers. Rows with stream order >= this value are plotted
#'   (inclusive).
#' @param bbox_buffer Numeric scalar or length-2 numeric vector. Buffer applied
#'   to the computed bounding box before plotting, in map units of the data CRS.
#'   If scalar, applied to both x and y directions. Default is `0` (no buffer).
#' @param style Named list of visual parameters. Only provided names override
#'   defaults. Use to avoid cluttering the function signature.
#' @param add_scale Logical. If `TRUE`, adds `ggspatial::annotation_scale()`.
#' @param add_north Logical. If `TRUE`, adds `ggspatial::annotation_north_arrow()`.
#'
#' @return A `ggplot` object with an additional attribute `"bbox"` containing
#' the unbuffered bounding box for all layers, and `"bbox_plot"` containing the
#' buffered bounding box actually used for plotting.
#'
#' @examples
#' \dontrun{
#' p <- plot_basin_map(
#'   basins_sf = basins_sf,
#'   subbasins_sf = subbasins_sf,
#'   rivers_sf = rivers_sf,
#'   points_sf = points_sf,
#'   min_stream_order = 3,
#'   bbox_buffer = 5000
#' )
#' print(p)
#' attr(p, "bbox")
#' attr(p, "bbox_plot")
#' }
#'
#' @importFrom rlang .data
#' @importFrom sf st_crs st_transform st_coordinates st_bbox
#' @importFrom dplyr filter
#' @import ggplot2
#' @export
plot_basin_map <- function(
    basins_sf = NULL,
    subbasins_sf = NULL,
    rivers_sf = NULL,
    points_sf = NULL,
    stream_order_col = "strord",
    label_col = "name",
    min_stream_order = 2,
    bbox_buffer = 5,
    style = list(),
    add_scale = FALSE,
    add_north = FALSE
) {
  # ---------------------------------------------------------------------------
  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (is.null(basins_sf) && is.null(subbasins_sf) && is.null(rivers_sf) && is.null(points_sf)) {
    stop(
      "At least one of `basins_sf`, `subbasins_sf`, `rivers_sf`, or `points_sf` must be provided.",
      call. = FALSE
    )
  }
  if (!is.null(basins_sf) && !inherits(basins_sf, "sf")) stop("`basins_sf` must be an sf object.", call. = FALSE)
  if (!is.null(subbasins_sf) && !inherits(subbasins_sf, "sf")) stop("`subbasins_sf` must be an sf object.", call. = FALSE)
  if (!is.null(rivers_sf) && !inherits(rivers_sf, "sf")) stop("`rivers_sf` must be an sf object.", call. = FALSE)
  if (!is.null(points_sf) && !inherits(points_sf, "sf")) stop("`points_sf` must be an sf object.", call. = FALSE)

  if (!is.character(stream_order_col) || length(stream_order_col) != 1L || is.na(stream_order_col)) {
    stop("`stream_order_col` must be a non-missing character scalar.", call. = FALSE)
  }
  if (!is.character(label_col) || length(label_col) != 1L || is.na(label_col)) {
    stop("`label_col` must be a non-missing character scalar.", call. = FALSE)
  }
  if (!is.numeric(min_stream_order) || length(min_stream_order) != 1L || is.na(min_stream_order)) {
    stop("`min_stream_order` must be a non-missing numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(bbox_buffer) || anyNA(bbox_buffer) || !(length(bbox_buffer) %in% c(1L, 2L))) {
    stop("`bbox_buffer` must be a non-missing numeric scalar or length-2 numeric vector.", call. = FALSE)
  }
  if (any(bbox_buffer < 0)) stop("`bbox_buffer` must be >= 0.", call. = FALSE)

  if (!is.list(style)) {
    stop("`style` must be a list.", call. = FALSE)
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
  } else if (!is.null(subbasins_sf)) {
    ref_crs <- sf::st_crs(subbasins_sf)
  } else if (!is.null(rivers_sf)) {
    ref_crs <- sf::st_crs(rivers_sf)
  } else {
    ref_crs <- sf::st_crs(points_sf)
  }

  if (is.null(ref_crs) || is.na(ref_crs)) {
    stop("Could not determine a valid CRS from inputs (missing/NA CRS).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # CRS harmonization
  # ---------------------------------------------------------------------------
  if (!is.null(basins_sf) && sf::st_crs(basins_sf) != ref_crs) basins_sf <- sf::st_transform(basins_sf, ref_crs)
  if (!is.null(subbasins_sf) && sf::st_crs(subbasins_sf) != ref_crs) subbasins_sf <- sf::st_transform(subbasins_sf, ref_crs)
  if (!is.null(rivers_sf) && sf::st_crs(rivers_sf) != ref_crs) rivers_sf <- sf::st_transform(rivers_sf, ref_crs)
  if (!is.null(points_sf) && sf::st_crs(points_sf) != ref_crs) points_sf <- sf::st_transform(points_sf, ref_crs)

  # ---------------------------------------------------------------------------
  # Compute combined bounding box for all layers (unbuffered)
  # ---------------------------------------------------------------------------
  bbox_list <- list()
  if (!is.null(basins_sf)) bbox_list$basins <- sf::st_bbox(basins_sf)
  if (!is.null(subbasins_sf)) bbox_list$subbasins <- sf::st_bbox(subbasins_sf)
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
  # Buffered bbox for plotting
  # ---------------------------------------------------------------------------
  buffer_xy <- if (length(bbox_buffer) == 1L) c(bbox_buffer, bbox_buffer) else bbox_buffer
  plot_bbox <- combined_bbox
  plot_bbox[["xmin"]] <- plot_bbox[["xmin"]] - buffer_xy[[1L]]
  plot_bbox[["xmax"]] <- plot_bbox[["xmax"]] + buffer_xy[[1L]]
  plot_bbox[["ymin"]] <- plot_bbox[["ymin"]] - buffer_xy[[2L]]
  plot_bbox[["ymax"]] <- plot_bbox[["ymax"]] + buffer_xy[[2L]]
  class(plot_bbox) <- "bbox"
  attr(plot_bbox, "crs") <- ref_crs

  # ---------------------------------------------------------------------------
  # Prepare point coordinates + label nudges (based on overall map extent)
  # ---------------------------------------------------------------------------
  points_df <- NULL
  nudge_x <- 0
  nudge_y <- 0

  if (!is.null(points_sf)) {
    if (!(label_col %in% names(points_sf))) {
      stop(sprintf("`points_sf` must contain a `%s` column for labelling.", label_col), call. = FALSE)
    }

    coords <- sf::st_coordinates(points_sf)
    points_df <- cbind(sf::st_drop_geometry(points_sf), coords)
    points_df[[".label"]] <- points_df[[label_col]]

    bbox_width  <- as.numeric(plot_bbox[["xmax"]] - plot_bbox[["xmin"]])
    bbox_height <- as.numeric(plot_bbox[["ymax"]] - plot_bbox[["ymin"]])

    NUDGE_X_FRAC <- 0.03
    NUDGE_Y_FRAC <- 0.04
    nudge_x <- if (is.finite(bbox_width) && bbox_width > 0) bbox_width * NUDGE_X_FRAC else 0
    nudge_y <- if (is.finite(bbox_height) && bbox_height > 0) bbox_height * NUDGE_Y_FRAC else 0
  }

  # ---------------------------------------------------------------------------
  # Plot styling defaults + override via `style`
  # ---------------------------------------------------------------------------
  style_defaults <- list(
    BASIN_FILL = "grey95",
    BASIN_BORDER = "grey70",
    BASIN_LINEWIDTH = 0.2,

    SUBBASIN_BORDER = "grey70",
    SUBBASIN_LINEWIDTH = 0.15,
    SUBBASIN_ALPHA = 1.0,

    RIVER_COLOR = "#4783A9",
    RIVER_ALPHA = 0.9,

    MAJOR_RIVER_WIDTH = 1.00,
    LARGE_RIVER_WIDTH = 0.60,
    MEDIUM_RIVER_WIDTH = 0.20,
    SMALL_RIVER_WIDTH = 0.04,

    LABEL_TEXT_SIZE = 3,
    LABEL_FILL = scales::alpha("white", 0.6),
    LABEL_COLOR = "black",

    POINT_SHAPE = 21,
    POINT_FILL = "black",
    POINT_BORDER = "white",
    POINT_SIZE = 3,
    POINT_STROKE = 1
  )

  style_names <- names(style)
  if (length(style) > 0 && (is.null(style_names) || any(!nzchar(style_names)))) {
    stop("`style` must be a named list when providing overrides.", call. = FALSE)
  }

  unknown_style <- if (length(style) == 0) character(0) else setdiff(style_names, names(style_defaults))
  if (length(unknown_style) > 0) {
    stop(sprintf("Unknown `style` names: %s", paste(unknown_style, collapse = ", ")), call. = FALSE)
  }

  style_use <- style_defaults
  if (length(style) > 0) style_use[style_names] <- style

  # ---------------------------------------------------------------------------
  # ggplot2 compatibility: linewidth aesthetic became reliable in 3.4.0+
  # ---------------------------------------------------------------------------
  use_linewidth <- TRUE
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    use_linewidth <- utils::packageVersion("ggplot2") >= "3.4.0"
  }

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot()

  if (!is.null(basins_sf)) {
    p <- p + geom_sf(
      data = basins_sf,
      fill = style_use$BASIN_FILL,
      color = style_use$BASIN_BORDER,
      linewidth = style_use$BASIN_LINEWIDTH
    )
  }

  if (!is.null(subbasins_sf)) {
    # Overlay only borders by default (no fill) to avoid masking basin fill/rivers.
    p <- p + geom_sf(
      data = subbasins_sf,
      fill = NA,
      color = style_use$SUBBASIN_BORDER,
      linewidth = style_use$SUBBASIN_LINEWIDTH,
      alpha = style_use$SUBBASIN_ALPHA
    )
  }

  if (!is.null(rivers_sf)) {
    if (!(stream_order_col %in% names(rivers_sf))) {
      stop(sprintf("`rivers_sf` must contain a `%s` column for stream order.", stream_order_col), call. = FALSE)
    }
    if (!is.numeric(rivers_sf[[stream_order_col]])) {
      stop(sprintf("`rivers_sf$%s` must be numeric.", stream_order_col), call. = FALSE)
    }

    rivers_sf[[".stream_order"]] <- rivers_sf[[stream_order_col]]

    rivers_filt <- dplyr::filter(
      rivers_sf,
      .data$.stream_order >= min_stream_order
    )

    if (nrow(rivers_filt) == 0) {
      warning(
        sprintf(
          "No rivers remain after filtering with `min_stream_order = %s`. Consider lowering the threshold.",
          min_stream_order
        ),
        call. = FALSE
      )
    } else {
      rivers_filt <- dplyr::mutate(
        rivers_filt,
        .linewidth_class = factor(
          dplyr::case_when(
            .data$.stream_order >= 5 ~ "major",
            .data$.stream_order == 4 ~ "large",
            .data$.stream_order == 3 ~ "medium",
            TRUE ~ "small"
          ),
          levels = c("small", "medium", "large", "major")
        )
      )

      width_values <- c(
        major  = style_use$MAJOR_RIVER_WIDTH,
        large  = style_use$LARGE_RIVER_WIDTH,
        medium = style_use$MEDIUM_RIVER_WIDTH,
        small  = style_use$SMALL_RIVER_WIDTH
      )

      if (isTRUE(use_linewidth)) {
        p <- p +
          geom_sf(
            data = rivers_filt,
            aes(linewidth = .data$.linewidth_class),
            color = style_use$RIVER_COLOR,
            alpha = style_use$RIVER_ALPHA
          ) +
          scale_linewidth_manual(values = width_values, guide = "none")
      } else {
        p <- p +
          geom_sf(
            data = rivers_filt,
            aes(size = .data$.linewidth_class),
            color = style_use$RIVER_COLOR,
            alpha = style_use$RIVER_ALPHA
          ) +
          scale_size_manual(values = width_values, guide = "none")
      }
    }
  }

  if (!is.null(points_df)) {
    p <- p +
      geom_label(
        data = points_df,
        aes(x = .data$X, y = .data$Y, label = .data$.label),
        size = style_use$LABEL_TEXT_SIZE,
        linewidth = 0,
        fill = style_use$LABEL_FILL,
        color = style_use$LABEL_COLOR,
        label.r = grid::unit(0.1, "lines"),
        nudge_y = nudge_y,
        nudge_x = nudge_x
      ) +
      geom_point(
        data = points_df,
        aes(x = .data$X, y = .data$Y),
        shape = style_use$POINT_SHAPE,
        fill = style_use$POINT_FILL,
        color = style_use$POINT_BORDER,
        size = style_use$POINT_SIZE,
        stroke = style_use$POINT_STROKE
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
  # Final layout: apply buffered bbox so annotations/labels have room
  # ---------------------------------------------------------------------------
  p <- p +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    coord_sf(
      xlim = c(plot_bbox[["xmin"]], plot_bbox[["xmax"]]),
      ylim = c(plot_bbox[["ymin"]], plot_bbox[["ymax"]]),
      expand = FALSE
    )

  # Attach bounding boxes as attributes
  attr(p, "bbox") <- combined_bbox
  attr(p, "bbox_plot") <- plot_bbox

  p
}
