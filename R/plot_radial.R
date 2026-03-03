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
#' @param radial_breaks Numeric vector of radial axis breaks (must include 0).
#'   If `NULL`, computed via `pretty(c(0, max(value)))`.
#' @param radial_limits Numeric vector of length 2 giving `c(lower, upper)` for
#'   the radial axis. If `NULL` (default), computed as
#'   `c(0, max(radial_breaks) + step)` where step is the spacing between
#'   consecutive breaks. Category labels are placed one step inside the upper
#'   limit.
#' @param start_angle Numeric scalar. Angle in degrees where the first category
#'   is placed. Default 90 (top).
#' @param clockwise Logical scalar. If `TRUE`, categories proceed clockwise.
#'   Default `FALSE`.
#' @param show_grid_circles Logical scalar. If `TRUE`, draws concentric circles
#'   at positive `radial_breaks`. Default `TRUE`.
#' @param show_value_labels Logical scalar. If `TRUE`, adds numeric labels above
#'   bars (rounded to 1 decimal). Default `FALSE`.
#' @param legend_position Character scalar. One of `"right"` or `"bottom"`.
#'   Controls placement of the fill legend. Default `"right"`.
#'
#' @return A list of styling parameters.
#'
#' @examples
#' # Use defaults
#' radial_plot(df, category, value)
#'
#' # Customize specific styles
#' radial_plot(df, category, value,
#'   style = radial_style(bar_alpha = 0.6, category_label_size = 4, clockwise = TRUE)
#' )
#'
#' @export
radial_style <- function(
    bar_width = 0.75,
    bar_alpha = 0.85,
    bar_border_color = "grey40",
    category_label_size = 3,
    tick_label_size = 2.5,
    value_label_size = 2.5,
    spoke_color = "grey60",
    spoke_linetype = "dashed",
    grid_color = "grey90",
    strip_text_size = 11,
    radial_breaks = NULL,
    radial_limits = NULL,
    start_angle = 30,
    clockwise = FALSE,
    show_grid_circles = TRUE,
    show_value_labels = FALSE,
    legend_position = c("bottom", "right")
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
  if (!is.null(radial_limits)) {
    if (!is.numeric(radial_limits) || length(radial_limits) != 2L || anyNA(radial_limits)) {
      stop("`radial_limits` must be a numeric vector of length 2 with no NAs.", call. = FALSE)
    }
    if (radial_limits[1] >= radial_limits[2]) {
      stop("`radial_limits[1]` must be less than `radial_limits[2]`.", call. = FALSE)
    }
  }

  if (!is.numeric(start_angle) || length(start_angle) != 1L || is.na(start_angle)) {
    stop("`start_angle` must be a single non-NA numeric.", call. = FALSE)
  }
  if (!is.logical(clockwise) || length(clockwise) != 1L || is.na(clockwise)) {
    stop("`clockwise` must be a single non-NA logical.", call. = FALSE)
  }
  if (!is.logical(show_grid_circles) || length(show_grid_circles) != 1L || is.na(show_grid_circles)) {
    stop("`show_grid_circles` must be a single non-NA logical.", call. = FALSE)
  }
  if (!is.logical(show_value_labels) || length(show_value_labels) != 1L || is.na(show_value_labels)) {
    stop("`show_value_labels` must be a single non-NA logical.", call. = FALSE)
  }

  legend_position <- match.arg(legend_position)

  if (!is.null(radial_breaks)) {
    if (!is.numeric(radial_breaks) || length(radial_breaks) < 2L || anyNA(radial_breaks)) {
      stop("`radial_breaks` must be numeric, length >= 2, no NAs.", call. = FALSE)
    }
    if (!0 %in% radial_breaks) stop("`radial_breaks` must include 0.", call. = FALSE)
  }

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
      strip_text_size = strip_text_size,
      radial_breaks = radial_breaks,
      radial_limits = radial_limits,
      start_angle = start_angle,
      clockwise = clockwise,
      show_grid_circles = show_grid_circles,
      show_value_labels = show_value_labels,
      legend_position = legend_position
    ),
    class = "radial_style"
  )
}


#' Create a faceted radial (spider/polar) bar plot
#'
#' @description
#' Generates a radial bar plot (polar-coordinate columns) for categorical groups
#' arranged around a circle, optionally faceted.
#'
#' @details
#' `category`, `value`, and `fill` are tidy-evaluated and may be unquoted column
#' names or a single character string naming a column in `data`. If `facet` is
#' `NULL`, all data are plotted in a single panel (internally labelled `"all"`).
#'
#' **Fill behavior** is inferred from the type of the `fill` column:
#' - *Discrete* (character or factor): bars are colored by level using
#'   `scale_fill_manual()`. `palette` should be a named vector mapping levels to
#'   colors, or an unnamed vector matched by position. A default green-to-red
#'   palette is used when `palette = NULL`.
#' - *Continuous* (numeric): bars are colored by a gradient using
#'   `scale_fill_gradient()` (2 colors) or `scale_fill_gradientn()` (3+ colors).
#'   `palette` should be a color vector of length >= 2. `limits` pins the gradient
#'   endpoints. A yellow-to-red default is used when `palette = NULL`.
#'
#' Missing values in `value` are dropped. Behavior is controlled by `na_handling`.
#' Spokes always extend to the outermost radial break. Layout and display options
#' are controlled via `radial_style()`.
#'
#' @param data A non-empty `data.frame`.
#' @param category Column for circular axis categories. Unquoted name or a single
#'   character string naming a column in `data`.
#' @param value Numeric column for bar length (radial magnitude). Must be numeric
#'   and not entirely `NA`.
#' @param facet Optional faceting column. If `NULL`, no faceting.
#'
#' @param fill Column to map to bar fill color. Defaults to `value` (continuous
#'   gradient by bar height). May point to a discrete (character/factor) or
#'   numeric column. Unquoted name or a single character string.
#' @param palette Colors for the fill scale.
#' @param limits Numeric vector of length 2 (`c(lo, hi)`). Only used when `fill`
#'   is continuous.
#' @param na_color Color used for `NA` fill values. Default `NA` (transparent).
#' @param facet_nrow Integer or `NULL`. Number of rows for `facet_wrap()`.
#' @param style A `radial_style()` object controlling aesthetics and layout.
#' @param na_handling Character scalar. One of `"warn"`, `"silent"`, `"error"`.
#'
#' @return A `ggplot2` plot object.
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
    fill  = NULL,
    palette  = NULL,
    limits   = NULL,
    na_color = NA,
    facet_nrow = NULL,
    style = radial_style(),
    na_handling = c("warn", "silent", "error")
) {

  na_handling <- match.arg(na_handling)

  if (!inherits(style, "radial_style")) style <- radial_style()
  sty <- style

  direction         <- if (sty$clockwise) -1 else 1
  radial_breaks     <- sty$radial_breaks
  radial_limits     <- sty$radial_limits
  start_angle       <- sty$start_angle
  show_grid_circles <- sty$show_grid_circles
  show_value_labels <- sty$show_value_labels
  legend_position   <- sty$legend_position

  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }

  # Legend bar geometry (pt)
  if (legend_position == "bottom") {
    legend_barwidth_pt  <- 150
    legend_barheight_pt <- 10
  } else {
    legend_barwidth_pt  <- 10
    legend_barheight_pt <- 150
  }

  # ---------------------------------------------------------------------------
  # Capture fill quosure early (must happen before any nested enquo calls)
  # ---------------------------------------------------------------------------
  fill_quo     <- rlang::enquo(fill)
  fill_is_null <- rlang::quo_is_null(fill_quo)

  # ---------------------------------------------------------------------------
  # Column resolution
  # ---------------------------------------------------------------------------
  .as_col <- function(x, data) {
    x_quo <- rlang::enquo(x)
    if (rlang::quo_is_call(x_quo) || rlang::quo_is_symbol(x_quo)) return(x_quo)
    x_val <- rlang::eval_tidy(x_quo)
    if (is.character(x_val) && length(x_val) == 1L && x_val %in% names(data)) {
      return(rlang::sym(x_val))
    }
    x_quo
  }

  cat_col  <- .as_col({{ category }}, data)
  val_col  <- .as_col({{ value }},    data)

  # fill defaults to the value column
  fill_col       <- if (fill_is_null) val_col else .as_col({{ fill }}, data)
  fill_col_label <- rlang::as_label(fill_col)

  facet_quo <- rlang::enquo(facet)
  has_facet <- !rlang::quo_is_null(facet_quo)
  if (has_facet) facet_col <- .as_col({{ facet }}, data)

  # ---------------------------------------------------------------------------
  # Build internal data frame
  # ---------------------------------------------------------------------------
  if (has_facet) {
    df <- data |>
      dplyr::transmute(
        .cat      = as.character(!!cat_col),
        .val      = !!val_col,
        .facet    = as.character(!!facet_col),
        .fill_raw = !!fill_col
      )
  } else {
    df <- data |>
      dplyr::transmute(
        .cat      = as.character(!!cat_col),
        .val      = !!val_col,
        .facet    = "all",
        .fill_raw = !!fill_col
      )
  }

  if (!is.numeric(df$.val)) stop("`value` column must be numeric.", call. = FALSE)

  n_na <- sum(is.na(df$.val))
  if (all(is.na(df$.val))) stop("`value` contains only NA values.", call. = FALSE)

  if (n_na > 0L) {
    msg <- sprintf("`value` contains %d NA(s); these will not be plotted.", n_na)
    switch(na_handling,
           warn   = warning(msg, call. = FALSE),
           error  = stop(msg,    call. = FALSE),
           silent = NULL
    )
  }

  cat_levels   <- unique(df$.cat)
  facet_levels <- unique(df$.facet)

  df$.cat   <- factor(df$.cat,   levels = cat_levels)
  df$.facet <- factor(df$.facet, levels = facet_levels)

  # ---------------------------------------------------------------------------
  # Infer fill type from the mapped column
  # ---------------------------------------------------------------------------
  fill_raw <- df$.fill_raw

  fill_is_discrete   <- is.character(fill_raw) || is.factor(fill_raw)
  fill_is_continuous <- is.numeric(fill_raw)

  if (!fill_is_discrete && !fill_is_continuous) {
    stop(
      sprintf("`fill` column ('%s') must be numeric, character, or factor.", fill_col_label),
      call. = FALSE
    )
  }

  # Legend title is always the fill column name
  legend_title <- fill_col_label

  # ---------------------------------------------------------------------------
  # Fill processing
  # ---------------------------------------------------------------------------
  if (fill_is_discrete) {

    # Determine level order: honour factor levels, else appearance order
    if (is.factor(fill_raw)) {
      fill_levels_ord <- levels(fill_raw)
    } else {
      fill_levels_ord <- as.character(unique(fill_raw[!is.na(fill_raw)]))
    }

    df$.fill <- factor(df$.fill_raw, levels = fill_levels_ord)

    # Resolve palette
    if (!is.null(palette)) {
      if (!is.character(palette)) {
        stop("`palette` must be a character vector of colors.", call. = FALSE)
      }
      if (!is.null(names(palette))) {
        missing_lvls <- setdiff(fill_levels_ord, names(palette))
        if (length(missing_lvls) > 0L) {
          stop(
            sprintf("`palette` is missing colors for level(s): %s",
                    paste(missing_lvls, collapse = ", ")
            ),
            call. = FALSE
          )
        }
        palette <- palette[fill_levels_ord]
      } else {
        if (length(palette) != length(fill_levels_ord)) {
          stop(
            sprintf("Unnamed `palette` must have length %d (one per fill level).",
                    length(fill_levels_ord)
            ),
            call. = FALSE
          )
        }
        names(palette) <- fill_levels_ord
      }
    } else {
      default_pal <- c("#2E7D32", "#F0E442", "#E69F00", "#C00000")
      n <- length(fill_levels_ord)
      pal <- if (n <= length(default_pal)) {
        default_pal[seq_len(n)]
      } else {
        grDevices::colorRampPalette(c(default_pal[1L], default_pal[length(default_pal)]))(n)
      }
      names(pal) <- fill_levels_ord
      palette <- pal
    }

    fill_na_val <- factor(NA_character_, levels = fill_levels_ord)

  } else {

    df$.fill <- df$.fill_raw

    if (!is.null(palette)) {
      if (!is.character(palette) || length(palette) < 2L) {
        stop(
          "For continuous fill, `palette` must be a character vector of length >= 2.",
          call. = FALSE
        )
      }
    } else {
      palette <- c("yellow", "red")
    }
    palette <- unname(palette)

    if (!is.null(limits)) {
      if (!is.numeric(limits) || length(limits) != 2L ||
          anyNA(limits) || !all(is.finite(limits))) {
        stop("`limits` must be a numeric vector of length 2 with finite values.", call. = FALSE)
      }
      if (limits[1L] >= limits[2L]) {
        stop("`limits[1]` must be less than `limits[2]`.", call. = FALSE)
      }
      lims <- limits
    } else {
      lims <- range(df$.fill, na.rm = TRUE)
      if (!is.finite(lims[1L]) || !is.finite(lims[2L]) || lims[1L] == lims[2L]) {
        stop(
          sprintf(
            "Continuous fill on '%s' requires a finite, non-degenerate range. Supply `limits` to set manually.",
            fill_col_label
          ),
          call. = FALSE
        )
      }
    }

    fill_na_val <- NA_real_
  }

  # ---------------------------------------------------------------------------
  # Radial scale
  # ---------------------------------------------------------------------------
  max_val <- max(df$.val, na.rm = TRUE)

  if (is.null(radial_breaks)) {
    radial_breaks <- pretty(c(0, max_val))
  } else {
    radial_breaks <- sort(unique(radial_breaks))
  }

  y_max <- max(radial_breaks)

  pos_breaks <- radial_breaks[radial_breaks > 0]
  radial_step <- if (length(pos_breaks) >= 2L) {
    min(diff(sort(pos_breaks)))
  } else {
    y_max
  }

  if (is.null(radial_limits)) {
    radial_limits <- c(0, y_max + radial_step)
  }

  y_outer <- radial_limits[2]
  label_y <- y_outer - radial_step / 3

  # ---------------------------------------------------------------------------
  # GAP sector (creates the opening gap in the polar plot)
  # ---------------------------------------------------------------------------
  gap_rows <- df |>
    dplyr::group_by(.data$.facet) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      .cat      = "GAP",
      .val      = 0,
      .fill_raw = fill_na_val,
      .fill     = fill_na_val
    )

  df_plot <- dplyr::bind_rows(df, gap_rows)
  df_plot$.cat <- factor(df_plot$.cat, levels = c(cat_levels, "GAP"))

  label_df <- data.frame(.cat = factor(cat_levels, levels = c(cat_levels, "GAP")))

  # Tick labels on the GAP axis
  tick_df <- expand.grid(
    .facet = facet_levels,
    .cat   = factor("GAP", levels = c(cat_levels, "GAP")),
    .y     = radial_breaks[radial_breaks > 0],
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Spokes: ALWAYS draw for every category (even when bar missing / .val is NA)
  # ---------------------------------------------------------------------------
  spoke_df <- expand.grid(
    .facet = facet_levels,
    .cat   = factor(cat_levels, levels = c(cat_levels, "GAP")),
    stringsAsFactors = FALSE
  )
  spoke_df$.yend <- y_max

  # Value-axis line at GAP position (so 5/10/15/20 read as an axis)
  axis_df <- data.frame(
    .facet = facet_levels,
    .cat   = factor("GAP", levels = c(cat_levels, "GAP")),
    .y     = 0,
    .yend  = y_max,
    stringsAsFactors = FALSE
  )

  if (show_value_labels) {
    value_label_df <- df_plot |>
      dplyr::filter(.data$.cat != "GAP", !is.na(.data$.val), .data$.val > 0)
  }

  # ---------------------------------------------------------------------------
  # Theme
  # ---------------------------------------------------------------------------
  base_theme <- theme_bw() +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = if (show_grid_circles) {
        element_line(color = sty$grid_color, linewidth = 0.3)
      } else {
        element_blank()
      },
      panel.spacing   = grid::unit(0.8, "lines"),
      legend.position = legend_position,
      legend.justification = "right",
      plot.margin     = margin(15, 15, 15, 15)
    )

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot(df_plot)

  # Spokes (always)
  p <- p + geom_segment(
    data = spoke_df,
    aes(x = .data$.cat, xend = .data$.cat, y = 0, yend = .data$.yend),
    linetype  = sty$spoke_linetype,
    color     = sty$spoke_color,
    linewidth = 0.4
  )

  # Value axis line (at GAP)
  # p <- p + geom_segment(
  #   data = axis_df,
  #   aes(x = .data$.cat, xend = .data$.cat, y = .data$.y, yend = .data$.yend),
  #   linetype  = sty$spoke_linetype,
  #   color     = sty$spoke_color,
  #   linewidth = 0.4
  # )

  bar_data <- dplyr::filter(df_plot, .data$.cat != "GAP", !is.na(.data$.val))

  p <- p + geom_col(
    data = bar_data,
    aes(x = .data$.cat, y = .data$.val, fill = .data$.fill),
    width       = sty$bar_width,
    color       = sty$bar_border_color,
    linewidth   = 0.3,
    show.legend = TRUE,
    alpha       = sty$bar_alpha
  )

  if (show_value_labels) {
    if (exists("value_label_df") && nrow(value_label_df) > 0L) {
      p <- p + geom_text(
        data = value_label_df,
        aes(x = .data$.cat, y = .data$.val, label = round(.data$.val, 1)),
        size     = sty$value_label_size,
        vjust    = -0.5,
        fontface = "bold"
      )
    }
  }

  p <- p + geom_label(
    data = label_df,
    aes(x = .data$.cat, y = label_y, label = .data$.cat),
    fill          = "white",
    linewidth     = 0,
    label.padding = grid::unit(0.15, "lines"),
    size          = sty$category_label_size,
    fontface      = "bold"
  )

  # Tick labels along the GAP column
  p <- p + geom_label(
    data = tick_df,
    aes(x = .data$.cat, y = .data$.y, label = .data$.y),
    size  = sty$tick_label_size,
    fill = "white",
    border.color = "white",
    color = "grey30"
  )

  p <- p +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      limits = c(radial_limits[1], y_outer),
      breaks = radial_breaks[radial_breaks > 0],
      expand = c(0, 0)
    )

  # ---------------------------------------------------------------------------
  # Fill scale
  # ---------------------------------------------------------------------------
  if (fill_is_discrete) {
    p <- p + scale_fill_manual(
      values       = palette,
      name         = legend_title,
      limits       = fill_levels_ord,
      drop         = FALSE,
      na.translate = FALSE,
      na.value     = na_color
    )
  } else {
    if (length(palette) == 2L) {
      p <- p + scale_fill_gradient(
        low      = palette[1L],
        high     = palette[2L],
        limits   = lims,
        oob      = scales::squish,
        name     = legend_title,
        na.value = na_color
      )
    } else {
      p <- p + scale_fill_gradientn(
        colors   = palette,
        limits   = lims,
        oob      = scales::squish,
        name     = legend_title,
        na.value = na_color
      )
    }
  }

  # Legend colorbar size (works for continuous; ignored by discrete guides)
  #p <- p + guides(
  #  fill = guide_colorbar(
  #    barwidth  = grid::unit(legend_barwidth_pt, "pt"),
  #    barheight = grid::unit(legend_barheight_pt, "pt")
  #  )
  #)

  p <- p + coord_polar(
    start     = (start_angle %% 360) * pi / 180,
    direction = direction,
    clip      = "off"
  )

  if (has_facet) {
    p <- p + facet_wrap(~ .facet, nrow = facet_nrow)
  }

  p + base_theme
}
