#' Styling options for radial_plot
#'
#' @description
#' Creates a list of visual styling parameters for `radial_plot()`.
#' All parameters have sensible defaults; only specify what you want to change.
#'
#' @param bar_width Width of bars (0.1-1). Default 0.75.
#' @param bar_alpha Transparency of bars (0-1). Default 0.85.
#' @param bar_border_color Border color for bars. Default "grey30".
#' @param category_label_size Size of category labels around rim. Default 3.
#' @param tick_label_size Size of radial tick labels. Default 2.5.
#' @param value_label_size Size of value labels on bars. Default 2.5.
#' @param spoke_color Color of spoke lines. Default "grey50".
#' @param spoke_linetype Linetype for spokes. Default "dashed".
#' @param grid_color Color of concentric grid circles. Default "grey80".
#' @param strip_text_size Size of facet strip text. Default 11.
#'
#' @return A list of styling parameters.
#'
#' @examples
#' # Use defaults
#' radial_plot(df, category, value, fill)
#'
#' # Customize specific styles
#' radial_plot(df, category, value, fill,
#'   style = radial_style(bar_alpha = 0.6, category_label_size = 4)
#' )
#'
#' @export
radial_style <- function(
    bar_width = 0.75,
    bar_alpha = 0.85,
    bar_border_color = "grey30",
    category_label_size = 3,
    tick_label_size = 2.5,
    value_label_size = 2.5,
    spoke_color = "grey50",
    spoke_linetype = "dashed",
    grid_color = "grey80",
    strip_text_size = 11
) {

  # Validation

  .check_range <- function(x, name, lo, hi) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < lo || x > hi) {
      stop(sprintf("`%s` must be a number between %s and %s.", name, lo, hi), call. = FALSE)
    }
  }

  .check_range(bar_width, "bar_width", 0.1, 1)
  .check_range(bar_alpha, "bar_alpha", 0, 1)
  .check_range(category_label_size, "category_label_size", 0.5, 10)
  .check_range(tick_label_size, "tick_label_size", 0.5, 10)
  .check_range(value_label_size, "value_label_size", 0.5, 10)
  .check_range(strip_text_size, "strip_text_size", 5, 20)

  structure(
    list(
      bar_width = bar_width,
      bar_alpha = bar_alpha,
      bar_border_color = bar_border_color,
      category_label_size = category_label_size,
      tick_label_size = tick_label_size,
      value_label_size = value_label_size,
      spoke_color = spoke_color,
      spoke_linetype = spoke_linetype,
      grid_color = grid_color,
      strip_text_size = strip_text_size
    ),
    class = "radial_style"
  )
}


#' Create a faceted radial (spider/polar) bar plot
#'
#' @description
#' Generates a radial bar plot for visualizing categorical groups arranged
#' around a circle with corresponding numeric values, optionally faceted.
#'
#' @param data A data frame.
#' @param category Column for circular axis categories (spokes).
#' @param value Numeric column for bar length (radial magnitude).
#' @param fill Column for fill color grouping. If NULL, uses a single default color.
#'   Can be discrete (factor/character) or continuous (numeric).
#' @param facet Column for faceting into panels. If NULL, no faceting.
#' @param fill_colors For discrete fill: named character vector of colors.
#'   For continuous fill: length-2 vector of low/high colors.
#'   If NULL, uses sensible defaults (viridis-inspired palette).
#' @param default_fill Single color used when `fill` is NULL. Default "#4575b4".
#' @param radial_breaks Numeric vector of radial axis breaks (must include 0).
#'   If NULL, computed automatically.
#' @param radial_expand Expansion factor beyond max break. Default 1.1.
#' @param start_angle Angle (degrees) where first category is placed. Default 90 (top).
#' @param clockwise Logical. If TRUE, categories proceed clockwise. Default FALSE.
#' @param show_spokes Logical. Show spoke lines? Default TRUE.
#' @param spokes_extend_to One of "max", "value", or "none". Default "max".
#' @param show_grid_circles Logical. Show concentric circles? Default TRUE.
#' @param show_value_labels Logical. Show value labels on bars? Default FALSE.
#' @param legend_title Title for fill legend. Default uses column name or NULL if no fill.
#' @param facet_nrow Number of rows for facet layout. NULL for auto.
#' @param theme_style One of "minimal", "clean", or "classic". Default "minimal".
#' @param style A `radial_style()` object for fine-grained visual control.
#' @param na_handling How to handle NA values: "warn", "silent", or "error".
#'
#' @return A ggplot object.
#'
#' @examples
#' df <- data.frame(
#'   site = rep(c("A", "B", "C", "D"), 3),
#'   scenario = rep(c("Low", "Mid", "High"), each = 4),
#'   value = runif(12, 5, 30),
#'   risk = sample(c("Low", "Medium", "High"), 12, replace = TRUE),
#'   score = runif(12, 0, 100)
#' )
#'
#' # No fill - single color
#' radial_plot(df, category = site, value = value)
#'
#' # Discrete fill with auto colors
#' radial_plot(df, category = site, value = value, fill = risk)
#'
#' # Continuous fill
#' radial_plot(df, site, value, fill = score)
#'
#' # With faceting and custom colors
#' radial_plot(df, site, value, risk,
#'   facet = scenario,
#'   fill_colors = c(Low = "green", Medium = "orange", High = "red")
#' )
#'
#' @import ggplot2
#' @importFrom rlang .data enquo as_label quo_is_null
#' @importFrom dplyr transmute mutate bind_rows filter group_by slice ungroup
#' @export
radial_plot <- function(
    data,
    category,
    value,
    facet = NULL,

    fill_mode = c("none", "binned", "continuous"),
    fill_levels = NULL,
    fill_colors = NULL,

    radial_breaks = NULL,
    radial_expand = 1.1,
    start_angle = 90,
    clockwise = FALSE,
    show_spokes = TRUE,
    spokes_extend_to = c("max", "value", "none"),
    show_grid_circles = TRUE,
    show_value_labels = FALSE,
    legend_title = NULL,
    facet_nrow = NULL,
    theme_style = c("minimal", "clean", "classic"),
    style = radial_style(),
    na_handling = c("warn", "silent", "error")
) {

  spokes_extend_to <- match.arg(spokes_extend_to)
  theme_style <- match.arg(theme_style)
  na_handling <- match.arg(na_handling)
  fill_mode <- match.arg(fill_mode)

  direction <- if (clockwise) -1 else 1

  if (!inherits(style, "radial_style")) style <- radial_style()
  sty <- style

  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }

  .validate_scalar <- function(x, name, type, valid_range = NULL) {
    if (!inherits(x, type) || length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` must be a single non-NA %s.", name, type), call. = FALSE)
    }
    if (!is.null(valid_range) && (x < valid_range[1] || x > valid_range[2])) {
      stop(sprintf("`%s` must be between %s and %s.", name, valid_range[1], valid_range[2]), call. = FALSE)
    }
  }

  .validate_scalar(start_angle, "start_angle", "numeric")
  .validate_scalar(radial_expand, "radial_expand", "numeric", c(1, 2))
  .validate_scalar(clockwise, "clockwise", "logical")
  .validate_scalar(show_spokes, "show_spokes", "logical")
  .validate_scalar(show_grid_circles, "show_grid_circles", "logical")
  .validate_scalar(show_value_labels, "show_value_labels", "logical")

  .as_col <- function(x, data) {
    x_quo <- rlang::enquo(x)
    if (rlang::quo_is_call(x_quo) || rlang::quo_is_symbol(x_quo)) return(x_quo)
    x_val <- rlang::eval_tidy(x_quo)
    if (is.character(x_val) && length(x_val) == 1L && x_val %in% names(data)) {
      return(rlang::sym(x_val))
    }
    x_quo
  }

  cat_col <- .as_col({{ category }}, data)
  val_col <- .as_col({{ value }}, data)

  facet_quo <- rlang::enquo(facet)
  has_facet <- !rlang::quo_is_null(facet_quo)
  if (has_facet) facet_col <- .as_col({{ facet }}, data)

  if (has_facet) {
    df <- data |>
      dplyr::transmute(
        .cat = as.character(!!cat_col),
        .val = !!val_col,
        .facet = as.character(!!facet_col)
      )
  } else {
    df <- data |>
      dplyr::transmute(
        .cat = as.character(!!cat_col),
        .val = !!val_col,
        .facet = "all"
      )
  }

  if (!is.numeric(df$.val)) stop("`value` column must be numeric.", call. = FALSE)

  n_na <- sum(is.na(df$.val))

  if (all(is.na(df$.val))) {
    stop("`value` contains only NA values.", call. = FALSE)
  }

  if (n_na > 0L) {
    msg <- sprintf("`value` contains %d NA(s); these will not be plotted.", n_na)
    switch(na_handling,
           warn = warning(msg, call. = FALSE),
           error = stop(msg, call. = FALSE),
           silent = NULL
    )
  }

  cat_levels <- unique(df$.cat)
  facet_levels <- unique(df$.facet)
  df$.cat <- factor(df$.cat, levels = cat_levels)
  df$.facet <- factor(df$.facet, levels = facet_levels)

  .yellow_red_palette <- function(n)
    grDevices::colorRampPalette(c("yellow", "red"))(n)

  if (is.null(legend_title)) {
    legend_title <- switch(
      fill_mode,
      none = NULL,
      binned = "Category",
      continuous = "Value"
    )
  }

  # ---------------------------------------------------------------------------
  # Fill handling
  # Key refactor:
  # - `fill_levels` are numeric breaks/stops only (names ignored).
  # - For binned fill, labels come from names(fill_colors) if provided; otherwise auto.
  # ---------------------------------------------------------------------------
  if (fill_mode == "none") {
    df$.fill <- factor("default")

  } else if (fill_mode == "binned") {

    if (is.null(fill_levels) || !is.numeric(fill_levels) || length(fill_levels) < 2L || anyNA(fill_levels)) {
      stop("For `fill_mode = \"binned\"`, `fill_levels` must be a numeric vector (length >= 2) with no NAs.", call. = FALSE)
    }

    brks <- sort(unique(as.numeric(fill_levels)))
    if (length(brks) < 2L) stop("`fill_levels` must contain at least 2 unique breakpoints.", call. = FALSE)

    n_bins <- length(brks) - 1L

    # Determine bin labels:
    # 1) If fill_colors is named and length matches bins -> use those names (preferred)
    # 2) Else -> auto labels [a,b)
    use_named_colors <- !is.null(fill_colors) &&
      is.character(fill_colors) &&
      !is.null(names(fill_colors)) &&
      length(fill_colors) == n_bins &&
      all(nzchar(names(fill_colors)))

    if (use_named_colors) {
      bin_labels <- names(fill_colors)
    } else {
      a <- brks[-length(brks)]
      b <- brks[-1]
      fmt <- function(x) ifelse(is.infinite(x), "Inf", format(x, trim = TRUE, scientific = FALSE))
      bin_labels <- paste0("[", fmt(a), ", ", fmt(b), ")")
    }

    # Cut into bins
    df$.fill <- cut(
      df$.val,
      breaks = brks,
      labels = bin_labels,
      include.lowest = TRUE,
      right = FALSE
    )
    df$.fill <- factor(df$.fill, levels = bin_labels)

    # Normalize fill_colors to a named vector aligned to bin_labels
    if (!is.null(fill_colors)) {
      if (!is.character(fill_colors)) stop("`fill_colors` must be a character vector of colors.", call. = FALSE)

      if (!is.null(names(fill_colors))) {
        # If named colors: names must match labels exactly
        if (!setequal(names(fill_colors), bin_labels) || length(fill_colors) != length(bin_labels)) {
          stop(
            "For binned fill, if `fill_colors` is named it must have exactly one color per bin, ",
            "and names must match the bin labels (i.e., the legend categories).",
            call. = FALSE
          )
        }
        fill_colors <- fill_colors[bin_labels]
      } else {
        # Unnamed: must match bins by position (allowed but less preferred)
        if (length(fill_colors) != length(bin_labels)) {
          stop(sprintf("Unnamed `fill_colors` must have length %d (one per bin).", length(bin_labels)), call. = FALSE)
        }
        names(fill_colors) <- bin_labels
      }
    } else {
      pal <- .yellow_red_palette(length(bin_labels))
      names(pal) <- bin_labels
      fill_colors <- pal
    }

  } else if (fill_mode == "continuous") {

    df$.fill <- df$.val

    if (is.null(fill_colors)) fill_colors <- c("yellow", "red")
    if (!is.character(fill_colors) || length(fill_colors) < 2L) {
      stop("For `fill_mode = \"continuous\"`, `fill_colors` must be a character vector with length >= 2.", call. = FALSE)
    }
    # Names are ignored for continuous gradients (avoid implying categoricals)
    fill_colors <- unname(fill_colors)

    if (!is.null(fill_levels)) {
      if (!is.numeric(fill_levels) || length(fill_levels) < 2L || anyNA(fill_levels)) {
        stop("For continuous fill, `fill_levels` must be numeric length >= 2 with no NAs.", call. = FALSE)
      }
      if (any(!is.finite(fill_levels))) {
        stop("For continuous fill, `fill_levels` must be finite (no Inf/-Inf).", call. = FALSE)
      }
      fill_levels <- sort(unique(as.numeric(fill_levels)))
      if (length(fill_levels) < 2L) stop("For continuous fill, `fill_levels` must contain at least 2 unique stop points.", call. = FALSE)
      lims <- range(fill_levels)
    } else {
      lims <- range(df$.fill, na.rm = TRUE)
      fill_levels <- lims
    }

    if (!is.finite(lims[1]) || !is.finite(lims[2]) || lims[1] == lims[2]) {
      stop("Continuous fill requires a finite, non-degenerate value range.", call. = FALSE)
    }

    # If exactly 2 colors but >2 stops, interpolate across stops
    if (length(fill_colors) == 2L && length(fill_levels) > 2L) {
      fill_colors <- grDevices::colorRampPalette(fill_colors)(length(fill_levels))
    }
    # If colors count doesn't match stops, rebuild ramp using first/last (predictable)
    if (length(fill_colors) != length(fill_levels)) {
      fill_colors <- grDevices::colorRampPalette(c(fill_colors[1], fill_colors[length(fill_colors)]))(length(fill_levels))
    }

    values01 <- (fill_levels - lims[1]) / (lims[2] - lims[1])
    values01 <- pmax(0, pmin(1, values01))
  }

  # ---------------------------------------------------------------------------
  # Radial scale
  # ---------------------------------------------------------------------------
  max_val <- max(df$.val, na.rm = TRUE)

  if (is.null(radial_breaks)) {
    radial_breaks <- pretty(c(0, max_val))
  } else {
    if (!is.numeric(radial_breaks) || length(radial_breaks) < 2L || anyNA(radial_breaks)) {
      stop("`radial_breaks` must be numeric, length >= 2, no NAs.", call. = FALSE)
    }
    if (!0 %in% radial_breaks) stop("`radial_breaks` must include 0.", call. = FALSE)
    radial_breaks <- sort(unique(radial_breaks))
  }

  y_max <- max(radial_breaks)
  y_outer <- y_max * radial_expand

  # GAP sector
  gap_rows <- df |>
    dplyr::group_by(.data$.facet) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      .cat = "GAP",
      .val = 0,
      .fill = if (fill_mode == "continuous") NA_real_ else NA_character_
    )

  df_plot <- dplyr::bind_rows(df, gap_rows)
  df_plot$.cat <- factor(df_plot$.cat, levels = c(cat_levels, "GAP"))

  label_df <- data.frame(.cat = factor(cat_levels, levels = c(cat_levels, "GAP")))

  tick_df <- expand.grid(
    .facet = facet_levels,
    .cat = factor("GAP", levels = c(cat_levels, "GAP")),
    .y = radial_breaks[radial_breaks > 0],
    stringsAsFactors = FALSE
  )

  if (show_spokes && spokes_extend_to != "none") {
    spoke_df <- df_plot |>
      dplyr::filter(.data$.cat != "GAP", !is.na(.data$.val)) |>
      dplyr::mutate(.yend = if (spokes_extend_to == "max") y_max else .data$.val)
  } else {
    spoke_df <- NULL
  }

  if (show_grid_circles) {
    circle_breaks <- radial_breaks[radial_breaks > 0 & radial_breaks <= y_max]
  }

  if (show_value_labels) {
    value_label_df <- df_plot |>
      dplyr::filter(.data$.cat != "GAP", !is.na(.data$.val), .data$.val > 0)
  }

  # Theme
  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.8, "lines"),
      legend.position = "right",
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )

  theme_mods <- switch(
    theme_style,
    minimal = ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = sty$strip_text_size)
    ),
    clean = ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = sty$strip_text_size),
      panel.background = ggplot2::element_rect(fill = "grey98", color = NA)
    ),
    classic = ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey85", color = "grey60"),
      strip.text = ggplot2::element_text(face = "bold", size = sty$strip_text_size)
    )
  )

  # Plot
  p <- ggplot2::ggplot(df_plot)

  if (show_grid_circles) {
    for (r in circle_breaks) {
      p <- p + ggplot2::geom_hline(yintercept = r, color = sty$grid_color, linewidth = 0.3)
    }
  }

  if (!is.null(spoke_df)) {
    p <- p + ggplot2::geom_segment(
      data = spoke_df,
      ggplot2::aes(x = .data$.cat, xend = .data$.cat, y = 0, yend = .data$.yend),
      linetype = sty$spoke_linetype,
      color = sty$spoke_color,
      linewidth = 0.4
    )
  }

  bar_data <- dplyr::filter(df_plot, .data$.cat != "GAP", !is.na(.data$.val))

  if (fill_mode == "none") {
    p <- p + ggplot2::geom_col(
      data = bar_data,
      ggplot2::aes(x = .data$.cat, y = .data$.val),
      fill = sty$default_fill,
      width = sty$bar_width,
      color = sty$bar_border_color,
      linewidth = 0.3,
      alpha = sty$bar_alpha
    )
  } else {
    p <- p + ggplot2::geom_col(
      data = bar_data,
      ggplot2::aes(x = .data$.cat, y = .data$.val, fill = .data$.fill),
      width = sty$bar_width,
      color = sty$bar_border_color,
      linewidth = 0.3,
      show.legend = TRUE,
      alpha = sty$bar_alpha
    )
  }

  if (show_value_labels) {
    if (exists("value_label_df") && nrow(value_label_df) > 0L) {
      p <- p + ggplot2::geom_text(
        data = value_label_df,
        ggplot2::aes(x = .data$.cat, y = .data$.val, label = round(.data$.val, 1)),
        size = sty$value_label_size,
        vjust = -0.5,
        fontface = "bold"
      )
    }
  }

  p <- p + ggplot2::geom_label(
    data = label_df,
    ggplot2::aes(x = .data$.cat, y = y_outer, label = .data$.cat),
    fill = "white",
    linewidth = 0,
    label.padding = grid::unit(0.15, "lines"),
    size = sty$category_label_size,
    fontface = "bold"
  )

  p <- p + ggplot2::geom_text(
    data = tick_df,
    ggplot2::aes(x = .data$.cat, y = .data$.y, label = .data$.y),
    size = sty$tick_label_size,
    color = "grey30"
  )

  p <- p +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, y_outer), expand = c(0, 0))

  if (fill_mode == "binned") {
    p <- p + ggplot2::scale_fill_manual(
      values = fill_colors,
      name = legend_title,
      limits = names(fill_colors),
      drop = FALSE,
      na.translate = FALSE
    )
  } else if (fill_mode == "continuous") {
    if (!is.null(fill_levels) && length(fill_levels) > 2L) {
      p <- p + ggplot2::scale_fill_gradientn(
        colors = fill_colors,
        values = values01,
        limits = lims,
        oob = scales::squish,
        name = legend_title,
        na.value = NA
      )
    } else {
      p <- p + ggplot2::scale_fill_gradient(
        low = fill_colors[1],
        high = fill_colors[length(fill_colors)],
        limits = lims,
        oob = scales::squish,
        name = legend_title,
        na.value = NA
      )
    }
  }

  p <- p + ggplot2::coord_polar(
    start = (start_angle %% 360) * pi / 180,
    direction = direction,
    clip = "off"
  )

  if (has_facet) {
    p <- p + ggplot2::facet_wrap(~ .facet, nrow = facet_nrow)
  }

  p + base_theme + theme_mods
}

