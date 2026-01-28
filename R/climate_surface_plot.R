# =============================================================================
# Climate Response Surface Visualization Functions
# =============================================================================
#
# Two main functions:
# - climate_surface_base(): Creates publication-quality contour surface plots
# - ggsave_smart(): Intelligent saving with journal/PPT sizing
#
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helper functions (package-level, unexported)
# -----------------------------------------------------------------------------

.assert_scalar_int_ge <- function(x, nm, min_val) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < min_val) {
    stop("'", nm, "' must be a single number >= ", min_val, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_scalar_num <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop("'", nm, "' must be a single finite numeric value.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_limits <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 2 || any(!is.finite(x)) || x[1] >= x[2]) {
    stop("'", nm, "' must be numeric length-2 with ", nm, "[1] < ", nm, "[2].", call. = FALSE)
  }
  invisible(TRUE)
}

.pretty_bracketed <- function(rng, n_target) {
  if (!is.numeric(rng) || length(rng) != 2 || any(!is.finite(rng)) || rng[1] >= rng[2]) {
    stop("Internal error: invalid 'rng' for contour breaks.", call. = FALSE)
  }
  if (!is.numeric(n_target) || length(n_target) != 1 || !is.finite(n_target) || n_target < 2) {
    stop("Internal error: invalid 'n_target' for contour breaks.", call. = FALSE)
  }

  lo0 <- rng[1]
  hi0 <- rng[2]
  span <- hi0 - lo0

  n_bins_target <- as.integer(max(1, round(n_target)))
  raw_step <- span / n_bins_target

  pow10 <- 10 ^ floor(log10(raw_step))
  frac <- raw_step / pow10

  nice_frac <- if (frac <= 1) {
    1
  } else if (frac <= 2) {
    2
  } else if (frac <= 2.5) {
    2.5
  } else if (frac <= 5) {
    5
  } else {
    10
  }
  step <- nice_frac * pow10

  lo <- floor(lo0 / step) * step
  hi <- ceiling(hi0 / step) * step

  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    stop("Internal error: failed to compute valid bracketed limits.", call. = FALSE)
  }

  br <- seq(lo, hi, by = step)
  dec <- max(0, -floor(log10(step)) + 2)
  br <- round(br, dec)
  br[1] <- lo
  br[length(br)] <- hi
  br
}

.make_diverging_bins <- function(n_bins, rng, thr, fail_dir,
                                 col_failure = "#df0000",
                                 col_safe = "#0033FF",
                                 col_mid = "#FFFFFF",
                                 col_fail_light = "#FEE5D9",
                                 col_safe_light = "#D6E3FF") {
  lo <- rng[1]
  hi <- rng[2]

  if (!is.finite(lo) || !is.finite(hi) || hi <= lo || n_bins <= 0) {
    return(rep(col_mid, max(1L, n_bins)))
  }

  thr_c <- min(max(thr, lo), hi)
  prop <- (thr_c - lo) / (hi - lo)
  prop <- max(0.05, min(0.95, prop))

  n_low <- max(1L, floor(n_bins * prop))
  n_high <- n_bins - n_low

  if (n_high < 1L) {
    n_high <- 1L
    n_low <- n_bins - 1L
  }
  if (n_low < 1L) {
    n_low <- 1L
    n_high <- n_bins - 1L
  }

  if (fail_dir == 1) {
    cols_low <- if (n_low == 1L) {
      col_mid
    } else {
      c(colorRampPalette(c(col_safe, col_safe_light))(n_low - 1L), col_mid)
    }
    cols_high <- if (n_high == 1L) {
      col_failure
    } else {
      colorRampPalette(c(col_fail_light, col_failure))(n_high)
    }
  } else {
    cols_low <- if (n_low == 1L) {
      col_mid
    } else {
      c(colorRampPalette(c(col_failure, col_fail_light))(n_low - 1L), col_mid)
    }
    cols_high <- if (n_high == 1L) {
      col_safe
    } else {
      colorRampPalette(c(col_safe_light, col_safe))(n_high)
    }
  }

  cols <- c(cols_low, cols_high)

  if (length(cols) < n_bins) {
    cols <- c(cols, rep(tail(cols, 1), n_bins - length(cols)))
  }
  if (length(cols) > n_bins) {
    cols <- cols[seq_len(n_bins)]
  }

  cols
}

.make_legend_labeller <- function(max_labels = 14L) {
  force(max_labels)
  function(breaks) {
    out <- rep("", length(breaks))
    finite <- is.finite(breaks)
    b <- breaks[finite]
    n <- length(b)

    if (n == 0L) return(out)

    if (n <= max_labels) {
      out[which(finite)] <- format(b, trim = TRUE, scientific = FALSE)
      return(out)
    }

    stride <- as.integer(ceiling(n / max_labels))
    idx_keep <- seq.int(1L, n, by = stride)

    if (idx_keep[length(idx_keep)] != n) {
      idx_keep <- c(idx_keep, n)
    }

    pos <- which(finite)[idx_keep]
    out[pos] <- format(b[idx_keep], trim = TRUE, scientific = FALSE)
    out
  }
}


# -----------------------------------------------------------------------------
# Main function: climate_surface_base
# -----------------------------------------------------------------------------

#' Create a climate response surface plot
#'
#' @param data A data frame containing the surface data.
#' @param x_var Character. Column name for x-axis variable.
#' @param y_var Character. Column name for y-axis variable.
#' @param z_var Character. Column name for response variable.
#' @param threshold Numeric. Critical threshold value. If NULL, inferred from baseline.
#' @param title Character. Plot title.
#' @param x_label Expression or character. X-axis label.
#' @param y_label Expression or character. Y-axis label.
#' @param x_suffix Character. Suffix for x-axis tick labels.
#' @param y_suffix Character. Suffix for y-axis tick labels.
#' @param failure_dir Integer. 1 = failure above threshold, -1 = failure below.
#' @param x_breaks Numeric vector. Custom x-axis breaks.
#' @param y_breaks Numeric vector. Custom y-axis breaks.
#' @param n_contours Integer. Target number of contour levels.
#' @param z_limits Numeric length-2. Optional limits for z-axis.
#' @param legend_barwidth_spec Numeric. Legend bar width (fraction if <=1, inches if >1).
#' @param legend_barheight_spec Numeric. Legend bar height (fraction if <=1, inches if >1).
#' @param text_size Numeric. Text size multiplier.
#' @param facet Logical. If TRUE, create faceted plot.
#' @param facet_by Character. Column name for faceting.
#' @param facet_levels Character vector. Subset and order of facet levels.
#'
#' @return A ggplot object with metadata attributes.
#' @export
climate_surface_base <- function(
    data = NULL,
    x_var = NULL,
    y_var = NULL,
    z_var = NULL,
    threshold = NULL,
    title = "Climate Response Surface",
    x_label = expression(Delta ~ "Precipitation"),
    y_label = expression(Delta ~ "Temperature"),
    x_suffix = "%",
    y_suffix = "\u00B0C",
    failure_dir = 1,
    x_breaks = NULL,
    y_breaks = NULL,
    n_contours = 15,
    z_limits = NULL,
    legend_barwidth_spec = 4.31,
    legend_barheight_spec = 0.25,
    text_size = 0.95,
    facet = FALSE,
    facet_by = NULL,
    facet_levels = NULL
) {


  # Constants
  BASELINE_TOL = 1e-6
  LEGEND_MAX_LABELS <- 14L

  COLOR_FAILURE <- "#df0000"
  COLOR_SAFE <- "#0033FF"
  COLOR_MID <- "#FFFFFF"
  COLOR_FAIL_LIGHT <- "#FEE5D9"
  COLOR_SAFE_LIGHT <- "#D6E3FF"

  # Input validation
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' must have at least one row.", call. = FALSE)
  }

  if (is.null(x_var) || is.null(y_var) || is.null(z_var)) {
    stop("'x_var', 'y_var', 'z_var' are required.", call. = FALSE)
  }
  if (!x_var %in% names(data)) {
    stop("Column '", x_var, "' not found in data.", call. = FALSE)
  }
  if (!y_var %in% names(data)) {
    stop("Column '", y_var, "' not found in data.", call. = FALSE)
  }
  if (!z_var %in% names(data)) {
    stop("Column '", z_var, "' not found in data.", call. = FALSE)
  }

  if (!failure_dir %in% c(1, -1)) {
    stop("'failure_dir' must be 1 or -1.", call. = FALSE)
  }

  .assert_scalar_int_ge(n_contours, "n_contours", 3)
  .assert_scalar_num(BASELINE_TOL, "BASELINE_TOL")

  if (!is.null(z_limits)) {
    .assert_limits(z_limits, "z_limits")
  }

  if (!is.numeric(legend_barwidth_spec) || length(legend_barwidth_spec) != 1 ||
      !is.finite(legend_barwidth_spec) || legend_barwidth_spec <= 0) {
    stop("'legend_barwidth_spec' must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(legend_barheight_spec) || length(legend_barheight_spec) != 1 ||
      !is.finite(legend_barheight_spec) || legend_barheight_spec <= 0) {
    stop("'legend_barheight_spec' must be a positive numeric scalar.", call. = FALSE)
  }

  if (isTRUE(facet)) {
    if (is.null(facet_by)) {
      stop("'facet_by' is required when facet=TRUE.", call. = FALSE)
    }
    if (!facet_by %in% names(data)) {
      stop("Column '", facet_by, "' not found in data.", call. = FALSE)
    }
  }

  # Subset columns
  keep_cols <- unique(c(x_var, y_var, z_var, if (isTRUE(facet)) facet_by))
  data <- data[, keep_cols, drop = FALSE]

  if (!is.numeric(data[[x_var]])) {
    stop("'", x_var, "' must be numeric.", call. = FALSE)
  }
  if (!is.numeric(data[[y_var]])) {
    stop("'", y_var, "' must be numeric.", call. = FALSE)
  }
  if (!is.numeric(data[[z_var]])) {
    stop("'", z_var, "' must be numeric.", call. = FALSE)
  }

  # Handle faceting
  if (isTRUE(facet) && !is.null(facet_levels)) {
    data <- dplyr::filter(data, .data[[facet_by]] %in% facet_levels)
    if (nrow(data) == 0L) {
      stop("No data remaining after filtering to 'facet_levels'.", call. = FALSE)
    }
    data[[facet_by]] <- factor(data[[facet_by]], levels = facet_levels)
  }

  # Axis breaks
  if (is.null(x_breaks)) {
    x_breaks <- sort(unique(data[[x_var]]))
  }
  if (is.null(y_breaks)) {
    y_breaks <- sort(unique(data[[y_var]]))
  }

  # Z range validation
  z_vals <- data[[z_var]]
  z_rng_data <- range(z_vals[is.finite(z_vals)], na.rm = TRUE)

  if (!is.finite(z_rng_data[1]) || !is.finite(z_rng_data[2])) {
    stop("z has no finite values.", call. = FALSE)
  }
  if (z_rng_data[1] == z_rng_data[2]) {
    stop("z has zero range; cannot create contours.", call. = FALSE)
  }

  z_rng <- if (is.null(z_limits)) z_rng_data else z_limits

  if (!is.null(z_limits)) {
    data[[z_var]] <- pmin(pmax(data[[z_var]], z_limits[1]), z_limits[2])
  }

  # Threshold inference
  if (is.null(threshold)) {
    baseline_idx <- which(
      abs(data[[x_var]]) <= BASELINE_TOL &
        abs(data[[y_var]]) <= BASELINE_TOL
    )
    if (length(baseline_idx) == 0) {
      stop(
        "Cannot infer threshold: no baseline point found (x~0, y~0). ",
        "Provide 'threshold' explicitly or adjust 'BASELINE_TOL'.",
        call. = FALSE
      )
    }
    threshold <- mean(data[[z_var]][baseline_idx], na.rm = TRUE)
    if (!is.finite(threshold)) {
      stop("Inferred threshold is not finite (baseline z values are all NA).", call. = FALSE)
    }
  }
  .assert_scalar_num(threshold, "threshold")

  if (!is.null(z_limits) && (threshold < z_limits[1] || threshold > z_limits[2])) {
    warning(
      "'threshold' (", signif(threshold, 4), ") is outside 'z_limits' [",
      z_limits[1], ", ", z_limits[2], "]; threshold line may not display.",
      call. = FALSE
    )
  }

  # Compute contour breaks and colors
  contour_breaks <- .pretty_bracketed(rng = z_rng, n_target = as.integer(n_contours))
  n_bins <- length(contour_breaks) - 1L

  bin_cols <- .make_diverging_bins(
    n_bins = n_bins,
    rng = z_rng,
    thr = threshold,
    fail_dir = failure_dir,
    col_failure = COLOR_FAILURE,
    col_safe = COLOR_SAFE,
    col_mid = COLOR_MID,
    col_fail_light = COLOR_FAIL_LIGHT,
    col_safe_light = COLOR_SAFE_LIGHT
  )

  bin_mid <- 0.5 * (contour_breaks[-1] + contour_breaks[-length(contour_breaks)])
  bin_vals <- scales::rescale(bin_mid, from = z_rng)

  eps <- 1e-9
  bin_vals <- cummax(bin_vals + seq_along(bin_vals) * eps)

  values_use <- c(0, pmin(pmax(bin_vals, 0), 1), 1)
  colors_use <- c(bin_cols[1], bin_cols, bin_cols[length(bin_cols)])

  ord <- order(values_use)
  values_use <- values_use[ord]
  colors_use <- colors_use[ord]

  keep <- !duplicated(values_use, fromLast = TRUE)
  values_use <- values_use[keep]
  colors_use <- colors_use[keep]

  values_use <- cummax(values_use + seq_along(values_use) * eps)
  values_use <- pmin(values_use, 1)

  # Theme and styling
  x_span <- diff(range(x_breaks))
  y_span <- diff(range(y_breaks))
  xy_ratio <- if (is.finite(x_span) && is.finite(y_span) && y_span > 0) {
    x_span / y_span
  } else {
    1
  }

  base_size <- 12 * text_size

  theme_surface <- theme_bw(base_size = base_size) +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      panel.grid         = element_blank(),
      plot.title         = element_text(size = base_size + 1.5, hjust = 0),
      axis.title         = element_text(size = base_size + 1),
      axis.text          = element_text(size = base_size - 1),
      axis.title.x       = element_text(margin = margin(t = 6)),
      axis.title.y       = element_text(margin = margin(r = 6)),
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.text        = element_text(size = base_size - 2.5),
      legend.box.spacing = grid::unit(0, "pt"),
      #legend.margin      = margin(0, 0, 0, 0),
      legend.box.margin  = margin(b = 6),
      plot.margin        = margin(15, 25, 15, 15, unit = "pt"),
      legend.ticks       = element_line(colour = "gray60", linewidth = 0.4),
      legend.ticks.length = grid::unit(2.5, "mm"),
      strip.background   = element_rect(fill = "grey90", color = "grey70"),
      strip.text         = element_text(size = base_size, face = "bold")
    )

  guide_fill <- guide_coloursteps(
    direction = "horizontal",
    barwidth  = grid::unit(legend_barwidth_spec, "in"),
    barheight = grid::unit(legend_barheight_spec, "in"),
    show.limits = TRUE,
    ticks = TRUE
  )

  # Build plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_contour_filled(
      aes(z = .data[[z_var]], fill = after_stat(level_mid)),
      breaks = contour_breaks
    ) +
    geom_contour(
      aes(z = .data[[z_var]]),
      breaks = threshold,
      color = "black",
      linewidth = 1
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = function(v) paste0(v, x_suffix),
      expand = c(0, 0),
      limits = range(x_breaks)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(v) paste0(v, y_suffix),
      expand = c(0, 0),
      limits = range(y_breaks)
    ) +
    scale_fill_stepsn(
      colors = colors_use,
      values = values_use,
      limits = range(contour_breaks),
      breaks = contour_breaks,
      labels = .make_legend_labeller(max_labels = LEGEND_MAX_LABELS),
      oob = scales::squish,
      guide = guide_fill
    ) +
    coord_fixed(ratio = xy_ratio, expand = FALSE) +
    labs(x = x_label, y = y_label, title = title, fill = "") +
    theme_surface

  # Faceting
  if (isTRUE(facet)) {
    p <- p +
      facet_wrap(vars(!!rlang::sym(facet_by))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
  }

  # Attach metadata
  attr(p, "legend_barwidth_spec")  <- legend_barwidth_spec
  attr(p, "legend_barheight_spec") <- legend_barheight_spec
  attr(p, "threshold")             <- threshold
  attr(p, "z_range")               <- z_rng
  attr(p, "contour_breaks")        <- contour_breaks
  attr(p, "x_var")  <- x_var
  attr(p, "y_var")  <- y_var

  p
}





#' Overlay GCM information on a climate response surface plot
#'
#' @description
#' Adds Global Climate Model (GCM) summary information—typically scenario–horizon
#' change signals—to an existing climate response surface created with
#' \code{\link{climate_surface_base}}.
#'
#' The function overlays:
#' \itemize{
#'   \item GCM points positioned in the same \code{(x, y)} climate-change space
#'         as the response surface (e.g. \eqn{\Delta}precipitation vs
#'         \eqn{\Delta}temperature),
#'   \item scenario-specific colors and horizon-specific point shapes,
#'   \item optional bivariate-normal confidence ellipses summarizing the GCM
#'         ensemble spread.
#' }
#'
#' The underlying response surface is **not modified**; this function only adds
#' new ggplot layers and legend entries.
#'
#' @details
#' \strong{Facet compatibility:} If the base surface plot is faceted
#' (via \code{facet_by} in \code{climate_surface_base}), the same facet variable
#' must be present in \code{gcm_data}. Otherwise, GCM points will be recycled
#' into all panels.
#'
#' \strong{Axis consistency:} When available, axis limits stored as attributes
#' on the base plot (\code{x_breaks}, \code{y_breaks}) are used to clip
#' \code{gcm_data} before plotting, ensuring GCM points and ellipses are confined
#' to the visible plotting domain.
#'
#' @param p
#' A ggplot object, typically returned by \code{\link{climate_surface_base}}.
#'
#' @param gcm_data
#' A data frame containing GCM summary points. At minimum, this must include
#' numeric columns for the x- and y-axes (e.g. precipitation and temperature
#' changes), and categorical columns defining scenarios and time horizons.
#'
#' @param x_var,y_var
#' Names of the x and y columns in \code{gcm_data}. If \code{NULL}, these are
#' inferred from attributes attached to \code{p} by
#' \code{climate_surface_base}.
#'
#' @param scenario_var
#' Name of the column in \code{gcm_data} identifying emission scenarios
#' (e.g. \code{"ssp126"}, \code{"ssp245"}).
#'
#' @param horizon_var
#' Name of the column in \code{gcm_data} identifying time horizons
#' (e.g. \code{"near"}, \code{"far"}).
#'
#' @param scenario_levels
#' Optional character vector specifying which scenarios to include and their
#' plotting order. If supplied, scenarios not listed are dropped.
#'
#' @param horizon_levels
#' Optional character vector specifying which horizons to include and their
#' plotting order.
#'
#' @param scenario_colors
#' Named character vector mapping scenario names to colors. Names must match the
#' values (or factor levels) in \code{scenario_var}.
#'
#' @param horizon_shapes
#' Named integer vector mapping horizon names to point shapes. Names must match
#' the values (or factor levels) in \code{horizon_var}.
#'
#' @param alpha
#' Transparency of GCM points (0–1). Useful for dense ensembles.
#'
#' @param size
#' Point size for GCM markers (ggplot2 units).
#'
#' @param stroke
#' Stroke width for GCM point outlines.
#'
#' @param show_legend
#' Logical. If \code{FALSE}, suppresses all GCM-related legends (scenario and
#' horizon), while leaving the surface legend untouched.
#'
#' @param ellipse_levels
#' Optional numeric vector of confidence levels (values strictly between 0 and 1)
#' at which to draw bivariate-normal ellipses using \code{stat_ellipse()}.
#' For example, \code{c(0.5, 0.9)}.
#'
#' @param ellipse_group
#' Controls grouping for ellipse estimation:
#' \itemize{
#'   \item \code{"none"}: one ellipse per level using all GCM points,
#'   \item \code{"scenario"}: separate ellipses per scenario,
#'   \item \code{"horizon"}: separate ellipses per horizon,
#'   \item \code{"scenario_horizon"}: separate ellipses for each
#'         scenario–horizon combination.
#' }
#'
#' @param ellipse_color
#' Line color for ellipses when \code{ellipse_group = "none"} or when a fixed
#' color is desired.
#'
#' @param ellipse_linetype
#' Line type for ellipses (e.g. \code{"dashed"}).
#'
#' @param ellipse_linewidth
#' Line width for ellipses.
#'
#' @return
#' A ggplot object consisting of the original climate response surface with
#' GCM points and optional ellipses overlaid.
#'
#' @seealso
#' \code{\link{climate_surface_base}}, \code{\link[ggplot2]{stat_ellipse}}
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang .data
#' @export
#' Overlay GCM information on a climate response surface plot
#'
#' Adds Global Climate Model (GCM) points (scenario × horizon) and optional
#' spread summaries (Gaussian or robust ellipses; optional KDE HDR contours)
#' on top of a climate response surface produced by climate_surface_base().
#'
#' Notes on methodology:
#' - ellipse_method = "norm": mean/covariance (Gaussian second-moment summary)
#' - ellipse_method = "robust": robust mean/covariance via MASS::cov.rob (MCD)
#' - KDE contours (kde_levels) are Highest Density Region (HDR) mass contours
#'   computed in standardized (x,y) space within each group, then mapped back.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang .data
#' @export
climate_surface_gcm_overlay <- function(
    p,
    gcm_data = NULL,
    x_var = NULL,
    y_var = NULL,
    scenario_var = "scenario",
    horizon_var  = "horizon",
    scenario_levels = NULL,
    horizon_levels  = NULL,
    scenario_colors = c(
      "ssp126" = "#003466",
      "ssp245" = "#f69320",
      "ssp370" = "#df0000",
      "ssp585" = "#980002"
    ),
    horizon_shapes = c("near" = 1, "far" = 4),
    alpha = 0.75,
    size = 1,
    stroke = 1,
    show_legend = TRUE,
    ellipse_levels = NULL,
    ellipse_group = c("none", "scenario", "horizon", "scenario_horizon"),
    ellipse_color = "gray30",
    ellipse_linetype = "dashed",
    ellipse_linewidth = 0.5,
    # --- NEW (backward-compatible) ---
    ellipse_method = c("norm", "robust"),
    kde_levels = NULL,
    kde_group = c("none", "scenario", "horizon", "scenario_horizon"),
    kde_color = "gray30",
    kde_linetype = "solid",
    kde_linewidth = 0.5,
    kde_n = 120L,
    kde_bw_adjust = 1.0
) {

  ellipse_group <- match.arg(ellipse_group)
  ellipse_method <- match.arg(ellipse_method)
  kde_group <- match.arg(kde_group)

  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------
  if (is.null(p) || !inherits(p, "ggplot")) stop("'p' must be a ggplot object.", call. = FALSE)
  if (is.null(gcm_data)) return(p)
  if (!is.data.frame(gcm_data)) stop("'gcm_data' must be a data.frame.", call. = FALSE)

  if (is.null(x_var)) x_var <- attr(p, "x_var")
  if (is.null(y_var)) y_var <- attr(p, "y_var")

  if (is.null(x_var) || is.null(y_var)) {
    stop("Provide 'x_var' and 'y_var' (or pass a plot from climate_surface_base() with these attributes).", call. = FALSE)
  }

  req <- c(x_var, y_var, scenario_var, horizon_var)
  miss <- setdiff(req, names(gcm_data))
  if (length(miss) > 0) stop("Missing columns in gcm_data: ", paste(miss, collapse = ", "), call. = FALSE)

  if (!is.numeric(gcm_data[[x_var]]) || !is.numeric(gcm_data[[y_var]])) {
    stop("gcm_data x/y columns must be numeric.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single number in [0,1].", call. = FALSE)
  }

  if (!is.null(ellipse_levels)) {
    if (!is.numeric(ellipse_levels) || any(!is.finite(ellipse_levels)) ||
        any(ellipse_levels <= 0 | ellipse_levels >= 1)) {
      stop("'ellipse_levels' must be numeric values strictly between 0 and 1.", call. = FALSE)
    }
  }

  if (!is.null(kde_levels)) {
    if (!is.numeric(kde_levels) || any(!is.finite(kde_levels)) ||
        any(kde_levels <= 0 | kde_levels >= 1)) {
      stop("'kde_levels' must be numeric values strictly between 0 and 1.", call. = FALSE)
    }
  }

  if (!is.numeric(kde_n) || length(kde_n) != 1 || !is.finite(kde_n) || kde_n < 40) {
    stop("'kde_n' must be a single integer >= 40.", call. = FALSE)
  }
  kde_n <- as.integer(kde_n)

  if (!is.numeric(kde_bw_adjust) || length(kde_bw_adjust) != 1 || !is.finite(kde_bw_adjust) || kde_bw_adjust <= 0) {
    stop("'kde_bw_adjust' must be a single number > 0.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helpers (local)
  # ---------------------------------------------------------------------------
  .get_group_id <- function(df, mode, scenario_var, horizon_var) {
    if (mode == "none") return(rep("all", nrow(df)))
    if (mode == "scenario") return(as.character(df[[scenario_var]]))
    if (mode == "horizon") return(as.character(df[[horizon_var]]))
    # scenario_horizon
    paste0(as.character(df[[scenario_var]]), " / ", as.character(df[[horizon_var]]))
  }

  .ellipse_df_from_mu_S <- function(mu, S, level, n = 181L) {
    # Ellipse for 2D covariance at given probability level (Gaussian interpretation)
    # radius^2 ~ ChiSq(df=2)
    r <- sqrt(stats::qchisq(level, df = 2))
    ee <- eigen(S, symmetric = TRUE)
    vals <- pmax(ee$values, 0)
    A <- ee$vectors %*% diag(sqrt(vals), 2, 2)
    t <- seq(0, 2 * pi, length.out = n)
    circle <- rbind(cos(t), sin(t))
    pts <- t(mu + A %*% circle)
    data.frame(x = pts[, 1], y = pts[, 2], stringsAsFactors = FALSE)
  }

  .robust_mu_S <- function(x, y) {
    # MASS is a recommended package (ships with R); avoids new dependencies.
    # Use MCD; returns robust center and covariance.
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for robust ellipses (MASS::cov.rob).", call. = FALSE)
    }
    m <- cbind(x, y)
    # Need at least 3 points; practical minimum ~5 for stability
    if (nrow(m) < 5L) return(NULL)
    fit <- MASS::cov.rob(m, method = "mcd")
    if (any(!is.finite(fit$center)) || any(!is.finite(fit$cov))) return(NULL)
    list(mu = as.numeric(fit$center), S = fit$cov)
  }

  .norm_mu_S <- function(x, y) {
    m <- cbind(x, y)
    if (nrow(m) < 3L) return(NULL)
    mu <- colMeans(m)
    S <- stats::cov(m)
    if (any(!is.finite(mu)) || any(!is.finite(S))) return(NULL)
    # covariance can be singular; eigen() handles, ellipse clamps negatives.
    list(mu = as.numeric(mu), S = S)
  }

  .kde_hdr_contours_df <- function(x, y, mass_levels, n = 120L, bw_adjust = 1.0) {
    # Compute HDR (highest density region) contours in *standardized* space.
    # Steps:
    # 1) standardize within group
    # 2) KDE on regular grid
    # 3) find density thresholds s.t. integral of densities above threshold = mass
    # 4) contourLines at those thresholds
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for KDE contours (MASS::kde2d).", call. = FALSE)
    }

    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 10L) return(NULL)

    sx <- stats::sd(x)
    sy <- stats::sd(y)
    if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) return(NULL)

    mx <- mean(x); my <- mean(y)
    xs <- (x - mx) / sx
    ys <- (y - my) / sy

    # Bandwidth in standardized units; nrd is ok as default + adjust
    # Use MASS::bandwidth.nrd (same as stats::bw.nrd0-ish) to keep dependencies minimal.
    bwx <- MASS::bandwidth.nrd(xs) * bw_adjust
    bwy <- MASS::bandwidth.nrd(ys) * bw_adjust
    if (!is.finite(bwx) || !is.finite(bwy) || bwx <= 0 || bwy <= 0) return(NULL)

    kd <- MASS::kde2d(xs, ys, n = n, h = c(bwx, bwy))

    # Approx integral weights
    dx <- diff(kd$x[1:2])
    dy <- diff(kd$y[1:2])
    if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) return(NULL)

    z <- kd$z
    if (!all(is.finite(z))) return(NULL)

    # HDR thresholds for each mass level
    # Sort densities, compute cumulative mass = sum(z)*dx*dy
    z_vec <- as.vector(z)
    ord <- order(z_vec, decreasing = TRUE)
    z_sorted <- z_vec[ord]
    cum_mass <- cumsum(z_sorted) * dx * dy
    total_mass <- sum(z_sorted) * dx * dy
    if (!is.finite(total_mass) || total_mass <= 0) return(NULL)
    cum_mass <- cum_mass / total_mass

    thresholds <- vapply(mass_levels, function(m) {
      i <- which(cum_mass >= m)[1]
      if (is.na(i)) min(z_sorted) else z_sorted[i]
    }, numeric(1))

    # Contour lines at those density thresholds
    cl_all <- list()
    for (i in seq_along(thresholds)) {
      lev <- thresholds[i]
      cl <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = lev)
      if (length(cl) == 0L) next
      for (j in seq_along(cl)) {
        df <- data.frame(
          x = cl[[j]]$x * sx + mx,
          y = cl[[j]]$y * sy + my,
          level = mass_levels[i],
          piece = paste0("p", i, "_", j),
          stringsAsFactors = FALSE
        )
        cl_all[[length(cl_all) + 1L]] <- df
      }
    }

    if (length(cl_all) == 0L) return(NULL)
    dplyr::bind_rows(cl_all)
  }

  # ---------------------------------------------------------------------------
  # Prepare data (drop NA, clip to base plot range when available)
  # ---------------------------------------------------------------------------
  g <- gcm_data %>%
    dplyr::filter(is.finite(.data[[x_var]]), is.finite(.data[[y_var]]))

  x_breaks <- attr(p, "x_breaks")
  y_breaks <- attr(p, "y_breaks")

  if (is.numeric(x_breaks) && length(x_breaks) >= 2) {
    g <- g %>% dplyr::filter(.data[[x_var]] >= min(x_breaks), .data[[x_var]] <= max(x_breaks))
  }
  if (is.numeric(y_breaks) && length(y_breaks) >= 2) {
    g <- g %>% dplyr::filter(.data[[y_var]] >= min(y_breaks), .data[[y_var]] <= max(y_breaks))
  }

  if (!is.null(scenario_levels)) {
    g <- g %>% dplyr::filter(.data[[scenario_var]] %in% scenario_levels)
    g[[scenario_var]] <- factor(g[[scenario_var]], levels = scenario_levels)
  }

  if (!is.null(horizon_levels)) {
    g <- g %>% dplyr::filter(.data[[horizon_var]] %in% horizon_levels)
    g[[horizon_var]] <- factor(g[[horizon_var]], levels = horizon_levels)
  }

  if (nrow(g) == 0) return(p)

  if (is.factor(g[[scenario_var]])) {
    need <- levels(g[[scenario_var]])
    miss_cols <- setdiff(need, names(scenario_colors))
    if (length(miss_cols) > 0) stop("scenario_colors missing: ", paste(miss_cols, collapse = ", "), call. = FALSE)
  }
  if (is.factor(g[[horizon_var]])) {
    need <- levels(g[[horizon_var]])
    miss_shapes <- setdiff(need, names(horizon_shapes))
    if (length(miss_shapes) > 0) stop("horizon_shapes missing: ", paste(miss_shapes, collapse = ", "), call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Points
  # ---------------------------------------------------------------------------
  p2 <- p +
    geom_point(
      data = g,
      mapping = aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        color = .data[[scenario_var]],
        shape = .data[[horizon_var]]
      ),
      alpha = alpha,
      size = size,
      stroke = stroke
    ) +
    scale_color_manual(values = scenario_colors) +
    scale_shape_manual(values = horizon_shapes) +
    labs(color = "Scenario", shape = "Horizon")

  # ---------------------------------------------------------------------------
  # Optional ellipses (norm via stat_ellipse OR robust via explicit ellipse geometry)
  # ---------------------------------------------------------------------------
  if (!is.null(ellipse_levels)) {

    if (ellipse_method == "norm") {

      if (ellipse_group == "none") {
        for (lv in ellipse_levels) {
          p2 <- p2 + stat_ellipse(
            data = g,
            mapping = aes(x = .data[[x_var]], y = .data[[y_var]]),
            level = lv,
            type = "norm",
            color = ellipse_color,
            linetype = ellipse_linetype,
            linewidth = ellipse_linewidth
          )
        }
      } else {
        grp <- .get_group_id(g, ellipse_group, scenario_var, horizon_var)
        g2 <- g
        g2[[".ellipse_group"]] <- grp

        for (lv in ellipse_levels) {
          p2 <- p2 + stat_ellipse(
            data = g2,
            mapping = aes(
              x = .data[[x_var]],
              y = .data[[y_var]],
              group = .data[[".ellipse_group"]]
            ),
            level = lv,
            type = "norm",
            color = ellipse_color,
            linetype = ellipse_linetype,
            linewidth = ellipse_linewidth
          )
        }
      }

    } else {
      # Robust ellipses: robust mu/cov then analytic ellipse.
      grp <- .get_group_id(g, ellipse_group, scenario_var, horizon_var)
      g2 <- g
      g2[[".ellipse_group"]] <- grp

      ell_list <- list()
      ug <- unique(g2[[".ellipse_group"]])

      for (ggg in ug) {
        sub <- g2[g2[[".ellipse_group"]] == ggg, , drop = FALSE]
        fit <- .robust_mu_S(sub[[x_var]], sub[[y_var]])
        if (is.null(fit)) next

        for (lv in ellipse_levels) {
          df <- .ellipse_df_from_mu_S(fit$mu, fit$S, level = lv, n = 181L)
          names(df) <- c("x", "y")
          df[[".ellipse_group"]] <- ggg
          df[["level"]] <- lv
          df[["piece"]] <- paste0("ell_", ggg, "_", format(lv))
          ell_list[[length(ell_list) + 1L]] <- df
        }
      }

      if (length(ell_list) > 0L) {
        ell_df <- dplyr::bind_rows(ell_list)
        p2 <- p2 + geom_path(
          data = ell_df,
          mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["piece"]]),
          inherit.aes = FALSE,
          color = ellipse_color,
          linetype = ellipse_linetype,
          linewidth = ellipse_linewidth
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Optional KDE HDR contours (nonparametric density mass regions)
  # ---------------------------------------------------------------------------
  if (!is.null(kde_levels)) {

    grp <- .get_group_id(g, kde_group, scenario_var, horizon_var)
    g2 <- g
    g2[[".kde_group"]] <- grp

    kde_list <- list()
    ug <- unique(g2[[".kde_group"]])

    for (ggg in ug) {
      sub <- g2[g2[[".kde_group"]] == ggg, , drop = FALSE]

      d <- .kde_hdr_contours_df(
        x = sub[[x_var]],
        y = sub[[y_var]],
        mass_levels = sort(unique(kde_levels)),
        n = kde_n,
        bw_adjust = kde_bw_adjust
      )
      if (is.null(d) || nrow(d) == 0L) next
      d[[".kde_group"]] <- ggg
      kde_list[[length(kde_list) + 1L]] <- d
    }

    if (length(kde_list) > 0L) {
      kde_df <- dplyr::bind_rows(kde_list)

      # one path per contour piece
      p2 <- p2 + geom_path(
        data = kde_df,
        mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["piece"]]),
        inherit.aes = FALSE,
        color = kde_color,
        linetype = kde_linetype,
        linewidth = kde_linewidth
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Legend control (keep surface legend; toggle GCM legends only)
  # ---------------------------------------------------------------------------
  if (isTRUE(show_legend)) {
    p2 <- p2 + guides(
      color = guide_legend(
        title = "Scenario",
        order = 2,
        position = "right",
        direction = "vertical",
        theme = theme(legend.margin = margin(0, 0, 0, 10))
      ),
      shape = guide_legend(
        title = "Horizon",
        order = 3,
        position = "right",
        direction = "vertical",
        theme = theme(legend.margin = margin(0, 0, 0, 10))
      )
    )
  } else {
    p2 <- p2 + guides(color = "none", shape = "none")
  }

  p2
}





# -----------------------------------------------------------------------------
# Smart save function: ggsave_smart
# -----------------------------------------------------------------------------

#' Intelligent ggplot2 Save with Target-Aware Sizing
#'
#' @param filename Character. Output file path (must have extension).
#' @param plot A ggplot object.
#' @param target Character. Either "journal" or "ppt".
#' @param panel_w_in,panel_h_in Numeric. Target panel dimensions in inches.
#' @param gap_w_in,gap_h_in Numeric. Gap between panels in inches.
#' @param margin_left_in,margin_right_in,margin_top_in,margin_bottom_in Numeric.
#'   Figure margins in inches (can be 0).
#' @param journal_width_in Named numeric vector with single and double column widths.
#' @param journal_mode Character. "single" or "double" column mode.
#' @param ppt_width_in,ppt_height_in Numeric. PowerPoint slide dimensions.
#' @param dpi Numeric or NULL. Resolution for raster formats.
#' @param bg Character. Background color.
#' @param device Graphics device. NULL = auto-detect from extension.
#' @param verbose Logical. Print sizing information.
#'
#' @return Invisibly, a list containing plot and metadata.
#' @export
ggsave_smart <- function(
    filename,
    plot,
    target = c("journal", "ppt"),
    panel_w_in = 2.2,
    panel_h_in = 2.0,
    gap_w_in = 0.18,
    gap_h_in = 0.18,
    margin_left_in = 0.7,
    margin_right_in = 0.7,
    margin_top_in = 0.7,
    margin_bottom_in = 0.7,
    journal_width_in = c(single = 3.4, double = 6.7),
    journal_mode = c("double", "single"),
    ppt_width_in = 13.33,
    ppt_height_in = 7.5,
    dpi = NULL,
    bg = "white",
    device = NULL,
    verbose = FALSE
) {
  target <- match.arg(target)
  journal_mode <- match.arg(journal_mode)

  # Input validation
  if (!is.character(filename) || length(filename) != 1 || nchar(filename) == 0) {
    stop("'filename' must be a non-empty character scalar.", call. = FALSE)
  }

  ext <- tolower(tools::file_ext(filename))
  if (ext == "") {
    stop("'filename' must have a file extension (e.g., .pdf, .png).", call. = FALSE)
  }

  if (!inherits(plot, "ggplot")) {
    stop("'plot' must be a ggplot object.", call. = FALSE)
  }

  out_dir <- dirname(filename)
  if (nchar(out_dir) > 0 && out_dir != "." && !dir.exists(out_dir)) {
    ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      stop("Failed to create output directory: ", out_dir, call. = FALSE)
    }
  }

  .is_num1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)

  .assert_num1_ge0 <- function(x, nm) {
    if (!.is_num1(x) || x < 0) {
      stop("'", nm, "' must be a single number >= 0.", call. = FALSE)
    }
  }

  .assert_num1_gt0 <- function(x, nm) {
    if (!.is_num1(x) || x <= 0) {
      stop("'", nm, "' must be a single number > 0.", call. = FALSE)
    }
  }

  .assert_num1_gt0(panel_w_in, "panel_w_in")
  .assert_num1_gt0(panel_h_in, "panel_h_in")
  .assert_num1_ge0(gap_w_in, "gap_w_in")
  .assert_num1_ge0(gap_h_in, "gap_h_in")
  .assert_num1_ge0(margin_left_in, "margin_left_in")
  .assert_num1_ge0(margin_right_in, "margin_right_in")
  .assert_num1_ge0(margin_top_in, "margin_top_in")
  .assert_num1_ge0(margin_bottom_in, "margin_bottom_in")

  if (!is.numeric(journal_width_in) || length(journal_width_in) != 2 ||
      any(!is.finite(journal_width_in))) {
    stop("'journal_width_in' must be numeric length-2: c(single=..., double=...).", call. = FALSE)
  }

  # Infer facet layout
  .infer_layout <- function(p) {
    built <- ggplot_build(p)
    lay <- built$layout$layout
    n_panels <- if (!is.null(lay$PANEL)) length(unique(lay$PANEL)) else nrow(lay)
    n_panels <- max(1L, as.integer(n_panels))

    nrow <- 1L
    ncol <- 1L
    facet_params <- tryCatch(p$facet$params, error = function(e) NULL)

    if (!is.null(facet_params)) {
      if (!is.null(facet_params$nrow) && is.numeric(facet_params$nrow) && facet_params$nrow > 0) {
        nrow <- as.integer(facet_params$nrow)
      }
      if (!is.null(facet_params$ncol) && is.numeric(facet_params$ncol) && facet_params$ncol > 0) {
        ncol <- as.integer(facet_params$ncol)
      }
    }

    if (n_panels > 1L) {
      if (nrow > 1L && ncol == 1L) {
        ncol <- as.integer(ceiling(n_panels / nrow))
      }
      if (ncol > 1L && nrow == 1L) {
        nrow <- as.integer(ceiling(n_panels / ncol))
      }
      if (nrow == 1L && ncol == 1L) {
        ncol <- as.integer(ceiling(sqrt(n_panels)))
        nrow <- as.integer(ceiling(n_panels / ncol))
      }
    }

    list(nrow = nrow, ncol = ncol, n_panels = n_panels)
  }

  layout <- .infer_layout(plot)
  nrow <- layout$nrow
  ncol <- layout$ncol
  n_panels <- layout$n_panels

  # Compute figure size
  .compute_size <- function(nr, nc, pw, ph) {
    width <- margin_left_in + margin_right_in + nc * pw + max(0, nc - 1L) * gap_w_in
    height <- margin_top_in + margin_bottom_in + nr * ph + max(0, nr - 1L) * gap_h_in
    list(width = width, height = height)
  }

  size0 <- .compute_size(nrow, ncol, panel_w_in, panel_h_in)
  width_in <- size0$width
  height_in <- size0$height

  # Target-specific adjustments
  if (is.null(device)) {
    device <- switch(
      ext,
      pdf  = grDevices::cairo_pdf,
      png  = "png",
      tiff = "tiff",
      tif  = "tiff",
      svg  = "svg",
      eps  = "eps",
      jpg  = "jpeg",
      jpeg = "jpeg",
      NULL
    )
  }

  if (target == "journal") {
    if (is.null(dpi)) {
      dpi <- if (ext %in% c("png", "tif", "tiff", "jpg", "jpeg")) 300 else NA_real_
    }

    cap_w <- unname(journal_width_in[[journal_mode]])
    if (!.is_num1(cap_w) || cap_w <= 0) {
      stop("Invalid journal width cap for mode '", journal_mode, "'.", call. = FALSE)
    }

    if (n_panels == 1L) {
      scale_fac <- cap_w / width_in
      width_in <- cap_w
      height_in <- height_in * scale_fac
      if (verbose) {
        message("Single-panel journal: width = ", cap_w, " in; height scaled by ", signif(scale_fac, 3))
      }
    } else if (width_in > cap_w) {
      fixed_w <- margin_left_in + margin_right_in + max(0, ncol - 1L) * gap_w_in
      avail_w <- cap_w - fixed_w

      if (avail_w <= 0) {
        stop("Margins and gaps exceed journal width cap; reduce margins or gaps.", call. = FALSE)
      }

      scale_fac <- avail_w / (ncol * panel_w_in)
      panel_w_in2 <- panel_w_in * scale_fac
      panel_h_in2 <- panel_h_in * scale_fac

      size1 <- .compute_size(nrow, ncol, panel_w_in2, panel_h_in2)
      width_in <- size1$width
      height_in <- size1$height

      if (verbose) {
        message("Faceted journal: panels scaled by factor ", signif(scale_fac, 3))
      }
    }
  } else {
    if (is.null(dpi)) dpi <- 300

    scale_w <- ppt_width_in / width_in
    scale_h <- ppt_height_in / height_in
    scale_fac <- min(scale_w, scale_h)

    if (scale_fac < 1) {
      width_in <- width_in * scale_fac
      height_in <- height_in * scale_fac
      if (verbose) {
        message("PPT: scaled down by factor ", signif(scale_fac, 3), " to fit slide")
      }
    } else {
      width_in <- ppt_width_in
      height_in <- ppt_height_in
      if (verbose) {
        message("PPT: using full slide dimensions ", ppt_width_in, " x ", ppt_height_in, " in")
      }
    }
  }

  if (!.is_num1(width_in) || !.is_num1(height_in) || width_in <= 0 || height_in <= 0) {
    stop("Computed invalid output size.", call. = FALSE)
  }

  # Compute panel area and legend dimensions
  panel_area_w_in <- max(1e-6, width_in - margin_left_in - margin_right_in)
  panel_area_h_in <- max(1e-6, height_in - margin_top_in - margin_bottom_in)

  .spec_to_in <- function(spec, ref_in, default_frac) {
    if (!is.numeric(spec) || length(spec) != 1 || !is.finite(spec) || spec <= 0) {
      return(default_frac * ref_in)
    }
    if (spec <= 1) {
      return(spec * ref_in)
    }
    spec
  }

  bw_spec <- attr(plot, "legend_barwidth_spec", exact = TRUE)
  bh_spec <- attr(plot, "legend_barheight_spec", exact = TRUE)

  barwidth_in <- .spec_to_in(bw_spec, panel_area_w_in, default_frac = 0.85)
  barheight_in <- .spec_to_in(bh_spec, panel_area_h_in, default_frac = 0.04)

  barwidth_in <- max(1.0, min(panel_area_w_in * 0.95, barwidth_in))
  barheight_in <- max(0.10, min(0.50, barheight_in))

  plot <- plot + guides(
    fill = guide_coloursteps(
      direction = "horizontal",
      barwidth  = grid::unit(barwidth_in, "in"),
      barheight = grid::unit(barheight_in, "in"),
      show.limits = TRUE,
      ticks = TRUE
    )
  )

  # Save with cairo fallback for PDF
  .ggsave_do <- function(dev) {
    args <- list(
      filename = filename,
      plot = plot,
      width = width_in,
      height = height_in,
      units = "in",
      bg = bg,
      device = dev
    )
    if (is.finite(dpi)) {
      args$dpi <- dpi
    }
    do.call(ggsave, args)
    TRUE
  }

  ok <- FALSE
  err_msg <- NULL

  if (identical(ext, "pdf") && identical(device, grDevices::cairo_pdf)) {
    ok <- tryCatch(
      .ggsave_do(grDevices::cairo_pdf),
      error = function(e) {
        err_msg <<- conditionMessage(e)
        FALSE
      }
    )
    if (!ok) {
      if (verbose) {
        message("cairo_pdf failed (", err_msg, "); falling back to grDevices::pdf()")
      }
      ok <- tryCatch(
        .ggsave_do(grDevices::pdf),
        error = function(e) {
          err_msg <<- conditionMessage(e)
          FALSE
        }
      )
    }
  } else {
    ok <- tryCatch(
      .ggsave_do(device),
      error = function(e) {
        err_msg <<- conditionMessage(e)
        FALSE
      }
    )
  }

  if (!ok) {
    stop("ggsave failed: ", err_msg, call. = FALSE)
  }

  if (verbose) {
    message("Saved: ", filename, " (", width_in, " x ", height_in, " in)")
  }

  invisible(list(
    plot = plot,
    filename = filename,
    target = target,
    nrow = nrow,
    ncol = ncol,
    n_panels = n_panels,
    width_in = width_in,
    height_in = height_in,
    dpi = dpi,
    device = device,
    panel_area_w_in = panel_area_w_in,
    panel_area_h_in = panel_area_h_in,
    legend_barwidth_in = barwidth_in,
    legend_barheight_in = barheight_in
  ))
}
