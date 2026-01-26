#' Generate a Climate Response Surface Base Plot
#'
#' @description
#' Creates the base **climate response surface** (stress-test grid) as a filled contour plot.
#' This function plots only the surface (impacts as a function of two climate drivers) and
#' does **not** overlay GCM points; use your overlay helper for that.
#'
#' The surface is rendered with discrete filled contour bins (via
#' \code{\link[ggplot2:geom_contour_filled]{geom_contour_filled()}}) and a diverging palette
#' centered on a user-defined or inferred \code{threshold}. A bold contour line is drawn at
#' the threshold to delineate the transition boundary.
#'
#' The legend and surface share the same \code{z_limits} range (if provided). Values outside
#' \code{z_limits} are clipped to the range before plotting, and the legend is squished to
#' the same limits.
#'
#' @param data A data.frame containing the stress-test grid.
#'   Must include numeric columns for \code{x_var}, \code{y_var}, and \code{z_var}. If
#'   \code{facet=TRUE}, must also include \code{facet_by}.
#' @param x_var,y_var,z_var Character scalars. Column names in \code{data} for the x-axis
#'   driver, y-axis driver, and response (impact) variable, respectively.
#' @param threshold Numeric scalar. Response threshold used to:
#'   (i) draw the boundary contour line, and (ii) place the single pure-white step in the
#'   diverging palette. If \code{NULL}, the function infers the threshold from the baseline
#'   point where \code{x_var \eqn{\approx} 0} and \code{y_var \eqn{\approx} 0}. If no baseline
#'   is present, \code{threshold} must be provided.
#' @param title Plot title.
#' @param x_label,y_label Axis titles. Can be strings or \code{\link[base:plotmath]{plotmath}}
#'   expressions (defaults are suitable for delta precipitation/temperature).
#' @param x_suffix,y_suffix Character scalars appended to x/y tick labels (e.g. \code{"%"},
#'   \code{"\u00B0C"}).
#' @param failure_dir Integer, either \code{1} or \code{-1}. Controls which side of the
#'   threshold is shown as "failure" (red) vs "safe" (blue):
#'   \itemize{
#'     \item \code{failure_dir = 1}: higher \code{z_var} relative to \code{threshold} tends
#'       toward "safe" (blue) and lower tends toward "failure" (red).
#'     \item \code{failure_dir = -1}: the mapping is reversed.
#'   }
#'   Use this to align color semantics with your impact definition (e.g., higher values can
#'   be worse for some metrics).
#' @param x_breaks,y_breaks Optional numeric vectors of tick/level values for x and y axes.
#'   If \code{NULL}, uses sorted unique values from \code{data[[x_var]]} and \code{data[[y_var]]}.
#'   Axis limits are set to \code{range(x_breaks)} and \code{range(y_breaks)}.
#' @param n_contours Integer \eqn{\ge 3}. Target number of contour bins (filled steps).
#'   Internally used to generate a set of bracketed "pretty" breaks; the endpoints always
#'   match the plotted \code{z_limits} range.
#' @param z_limits Optional numeric length-2 vector \code{c(zmin, zmax)} with \code{zmin < zmax}.
#'   If provided, controls both:
#'   \itemize{
#'     \item the surface range (data are clipped to these limits), and
#'     \item the legend range (same limits, with out-of-bounds values squished).
#'   }
#'   If \code{NULL}, uses the finite range of \code{data[[z_var]]}.
#' @param panel_size_in Numeric scalar. Intended width (inches) used to size the legend bar.
#'   (The plot panel sizing returned in \code{$width} and \code{$height} is handled separately.)
#' @param legend_barheight_in Numeric scalar. Legend bar height in inches.
#' @param text_size Numeric scalar scaling factor applied to theme text sizes.
#' @param facet Logical. If \code{TRUE}, facet the surface using \code{facet_wrap()} by
#'   \code{facet_by}.
#' @param facet_by Character scalar. Column name used for faceting when \code{facet=TRUE}.
#' @param facet_levels Optional character vector specifying facet level order and/or subset.
#'   When provided, \code{data} are filtered to these levels and coerced to a factor with
#'   \code{levels = facet_levels}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{p}}{A \code{ggplot} object (filled contour surface).}
#'   \item{\code{width}}{Recommended output width in inches (numeric scalar).}
#'   \item{\code{height}}{Recommended output height in inches (numeric scalar).}
#' }
#'
#' The returned ggplot object \code{p} includes metadata attributes used by overlay/helper
#' functions:
#' \itemize{
#'   \item \code{threshold}
#'   \item \code{x_breaks}, \code{y_breaks}
#'   \item \code{z_limits}
#'   \item \code{contour_breaks}
#'   \item \code{x_var}, \code{y_var}, \code{z_var}
#' }
#'
#' @details
#' \strong{Threshold inference.} If \code{threshold=NULL}, the function searches for baseline
#' points where \code{x_var} and \code{y_var} are approximately zero (tolerance \eqn{1e-10})
#' and uses the mean \code{z_var} at those points. If none exist, an error is raised.
#'
#' \strong{Color scaling.} The diverging palette is constructed so that there is exactly one
#' pure-white step at the threshold boundary. The mapping flips depending on \code{failure_dir}.
#'
#' \strong{Legend labels.} The legend uses contour breaks; labels are thinned if necessary but
#' always include the minimum and maximum break, and a label at (or immediately below) the
#' threshold break.
#'
#' @examples
#' \dontrun{
#' grid <- expand.grid(
#'   dP = seq(-30, 30, by = 10),
#'   dT = seq(-2, 4, by = 1)
#' )
#' grid$impact <- 100 + 2 * grid$dT - 0.4 * grid$dP
#'
#' out <- climate_surface_base(
#'   data = grid,
#'   x_var = "dP",
#'   y_var = "dT",
#'   z_var = "impact",
#'   threshold = 100,
#'   n_contours = 21,
#'   z_limits = c(80, 120)
#' )
#'
#' out$p
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang .data sym
#' @importFrom scales rescale squish
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
    n_contours = 25,
    z_limits = NULL,
    panel_size_in = 6.0,
    legend_barheight_in = 0.25,
    text_size = 0.7,
    facet = FALSE,
    facet_by = NULL,
    facet_levels = NULL
) {

  # ---------------------------------------------------------------------------
  # Constants
  # ---------------------------------------------------------------------------
  BASELINE_TOL <- 1e-10
  Z_EXPAND     <- 1e-6
  LEGEND_MAX_LABELS <- 14L

  COLOR_FAILURE    <- "#df0000"
  COLOR_SAFE       <- "#0033FF"
  COLOR_MID        <- "#FFFFFF"
  COLOR_FAIL_LIGHT <- "#FEE5D9"
  COLOR_SAFE_LIGHT <- "#D6E3FF"

  panel_in = 5.5
  left_right_in <- 1.5
  bottom_in <- 1.0
  top_in <- 1.0

  # ---------------------------------------------------------------------------
  # Internal helpers
  # ---------------------------------------------------------------------------
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

  # Exact-number contour breaks (n_bins == n_contours)
  # Expanded, nice contour breaks
  .pretty_bracketed <- function(rng, n_target) {

    if (!is.numeric(rng) || length(rng) != 2 ||
        any(!is.finite(rng)) || rng[1] >= rng[2]) {
      stop("Internal error: invalid 'rng' for contour breaks.", call. = FALSE)
    }
    if (!is.numeric(n_target) || length(n_target) != 1 ||
        !is.finite(n_target) || n_target < 2) {
      stop("Internal error: invalid 'n_target' for contour breaks.", call. = FALSE)
    }

    lo0 <- rng[1]
    hi0 <- rng[2]

    # --- choose a nice step from the span / n_target ---
    span <- hi0 - lo0
    raw_step <- span / n_target

    pow  <- 10 ^ floor(log10(raw_step))
    frac <- raw_step / pow
    step <- (if (frac <= 1) 1 else if (frac <= 2) 2 else if (frac <= 2.5) 2.5 else if (frac <= 5) 5 else 10) * pow

    # --- expand toa nice bounds ---
    lo <- floor(lo0 / step) * step
    hi <- ceiling(hi0 / step) * step

    # --- anchor to zero if the whole range is on one side ---
    if (lo0 >= 0) lo <- 0
    if (hi0 <= 0) hi <- 0

    # --- generate breaks, with floating-point-safe rounding ---
    br <- seq(lo, hi, by = step)

    # round to a sensible number of decimals based on step (avoids 1090.0900000001)
    dec <- max(0, -floor(log10(step)) + 2)
    br <- round(br, dec)

    br[1] <- lo
    br[length(br)] <- hi
    br
  }


  # Diverging palette with exactly one pure-white step at threshold boundary
  .make_diverging_bins_single_white <- function(n_bins, rng, thr, fail_dir) {
    lo <- rng[1]; hi <- rng[2]
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo || n_bins <= 0) {
      return(rep(COLOR_MID, max(1L, n_bins)))
    }

    thr_c <- min(max(thr, lo), hi)
    prop  <- (thr_c - lo) / (hi - lo)

    n_low  <- max(1L, floor(n_bins * prop))
    n_high <- n_bins - n_low
    if (n_high < 1L) { n_high <- 1L; n_low <- n_bins - 1L }
    if (n_low  < 1L) { n_low  <- 1L; n_high <- n_bins - 1L }

    if (fail_dir == 1) {
      cols_low  <- if (n_low  == 1L) COLOR_FAILURE else colorRampPalette(c(COLOR_FAILURE, COLOR_FAIL_LIGHT))(n_low)
      cols_high <- if (n_high == 1L) COLOR_MID else c(COLOR_MID, colorRampPalette(c(COLOR_SAFE_LIGHT, COLOR_SAFE))(n_high - 1L))
    } else {
      cols_low  <- if (n_low  == 1L) COLOR_SAFE else colorRampPalette(c(COLOR_SAFE, COLOR_SAFE_LIGHT))(n_low)
      cols_high <- if (n_high == 1L) COLOR_MID else c(COLOR_MID, colorRampPalette(c(COLOR_FAIL_LIGHT, COLOR_FAILURE))(n_high - 1L))
    }

    cols <- c(cols_low, cols_high)
    if (length(cols) < n_bins) cols <- c(cols, rep(tail(cols, 1), n_bins - length(cols)))
    if (length(cols) > n_bins) cols <- cols[seq_len(n_bins)]
    cols
  }

  # Labels:
  # - always show min, "white" (floor to previous break), max
  # - if more than max_labels, thin evenly but keep mandatory labels
  .make_legend_labeller <- function(max_labels = LEGEND_MAX_LABELS, white_value = NULL) {
    force(max_labels); force(white_value)

    function(breaks) {
      # Always return same length as 'breaks'
      out <- rep("", length(breaks))

      finite <- is.finite(breaks)
      b <- breaks[finite]
      n <- length(b)
      if (n == 0L) return(out)

      # Indices are within the finite subset 'b'
      idx_keep <- c(1L, n)

      # white tick: choose break at or immediately BEFORE white_value
      i_white <- integer(0)
      if (!is.null(white_value) && is.finite(white_value)) {
        cand <- which(b <= white_value)
        i_white <- if (length(cand) == 0L) 1L else max(cand)
        idx_keep <- unique(c(idx_keep, i_white))
      }

      # If small enough, label all finite breaks
      if (n <= max_labels) {
        out[which(finite)] <- format(b, trim = TRUE, scientific = FALSE)
        return(out)
      }

      # Add evenly spaced indices, then trim to <= max_labels keeping mandatory
      idx_even <- unique(round(seq(1, n, length.out = max(3L, max_labels))))
      idx_keep <- sort(unique(c(idx_keep, idx_even)))

      if (length(idx_keep) > max_labels) {
        mandatory <- sort(unique(c(1L, n, i_white)))
        cand2 <- setdiff(idx_keep, mandatory)
        keep_n_cand <- max_labels - length(mandatory)

        if (keep_n_cand <= 0) {
          idx_keep <- mandatory
        } else if (length(cand2) > keep_n_cand) {
          cand2 <- cand2[unique(round(seq(1, length(cand2), length.out = keep_n_cand)))]
          idx_keep <- sort(unique(c(mandatory, cand2)))
        } else {
          idx_keep <- sort(unique(c(mandatory, cand2)))
        }
      }

      # Write labels back into full-length output
      pos <- which(finite)[idx_keep]
      out[pos] <- format(b[idx_keep], trim = TRUE, scientific = FALSE)
      out
    }
  }


  .compute_facet_dims <- function(k,
                                  panel_w_in = 5.5,
                                  panel_h_in = 5.5,
                                  left_right_in = 1.5,
                                  top_in = 1.0,
                                  bottom_in = 1.0,
                                  gap_w_in = 0.25,
                                  gap_h_in = 0.25,
                                  ncol = NULL,
                                  nrow = NULL,
                                  max_cols = 3) {

    if (is.null(ncol) && is.null(nrow)) {
      ncol <- min(max_cols, ceiling(sqrt(k)))
      nrow <- ceiling(k / ncol)
    } else if (is.null(nrow)) {
      nrow <- ceiling(k / ncol)
    } else if (is.null(ncol)) {
      ncol <- ceiling(k / nrow)
    }

    width  <- left_right_in + ncol * panel_w_in + (ncol - 1) * gap_w_in
    height <- bottom_in + top_in + nrow * panel_h_in + (nrow - 1) * gap_h_in

    list(width = width, height = height, nrow = nrow, ncol = ncol)
  }


  # ---------------------------------------------------------------------------
  # Validation & preparation
  # ---------------------------------------------------------------------------
  if (is.null(data) || !is.data.frame(data)) stop("'data' must be a data.frame.", call. = FALSE)
  if (is.null(x_var) || is.null(y_var) || is.null(z_var)) stop("'x_var', 'y_var', 'z_var' are required.", call. = FALSE)
  if (!x_var %in% names(data)) stop("Column '", x_var, "' not found in data.", call. = FALSE)
  if (!y_var %in% names(data)) stop("Column '", y_var, "' not found in data.", call. = FALSE)
  if (!z_var %in% names(data)) stop("Column '", z_var, "' not found in data.", call. = FALSE)
  if (!failure_dir %in% c(1, -1)) stop("'failure_dir' must be 1 or -1.", call. = FALSE)
  .assert_scalar_int_ge(n_contours, "n_contours", 3)

  if (!is.null(z_limits)) .assert_limits(z_limits, "z_limits")

  if (isTRUE(facet)) {
    if (is.null(facet_by)) stop("'facet_by' is required when facet=TRUE.", call. = FALSE)
    if (!facet_by %in% names(data)) stop("Column '", facet_by, "' not found in data.", call. = FALSE)
  }

  keep_cols <- unique(c(x_var, y_var, z_var, if (isTRUE(facet)) facet_by))
  data <- data[, keep_cols, drop = FALSE]

  if (!is.numeric(data[[x_var]])) stop("'", x_var, "' must be numeric.", call. = FALSE)
  if (!is.numeric(data[[y_var]])) stop("'", y_var, "' must be numeric.", call. = FALSE)
  if (!is.numeric(data[[z_var]])) stop("'", z_var, "' must be numeric.", call. = FALSE)

  if (isTRUE(facet) && !is.null(facet_levels)) {
    data <- dplyr::filter(data, .data[[facet_by]] %in% facet_levels)
    data[[facet_by]] <- factor(data[[facet_by]], levels = facet_levels)
  }

  if (is.null(x_breaks)) x_breaks <- sort(unique(data[[x_var]]))
  if (is.null(y_breaks)) y_breaks <- sort(unique(data[[y_var]]))

  z_vals <- data[[z_var]]
  z_rng_data <- range(z_vals[is.finite(z_vals)], na.rm = TRUE)
  if (!is.finite(z_rng_data[1]) || !is.finite(z_rng_data[2])) stop("z has no finite values.", call. = FALSE)
  if (z_rng_data[1] == z_rng_data[2]) stop("z has zero range; cannot contour.", call. = FALSE)

  # One range for both surface and legend: use z_limits if supplied, else data range
  z_rng <- if (is.null(z_limits)) z_rng_data else z_limits

  # Clip to range so surface respects z_rng
  if (!is.null(z_limits)) {
    data[[z_var]] <- pmin(pmax(data[[z_var]], z_limits[1]), z_limits[2])
  }

  # Threshold inference
  if (is.null(threshold)) {
    baseline_idx <- which(abs(data[[x_var]]) <= BASELINE_TOL & abs(data[[y_var]]) <= BASELINE_TOL)
    if (length(baseline_idx) == 0) {
      stop(
        "Cannot infer threshold: no baseline point (x≈0,y≈0). Provide 'threshold'.",
        call. = FALSE
      )
    }
    threshold <- mean(data[[z_var]][baseline_idx], na.rm = TRUE)
  }
  .assert_scalar_num(threshold, "threshold")

  # ---------------------------------------------------------------------------
  # Surface bins + colors
  # ---------------------------------------------------------------------------
  contour_breaks <- .pretty_bracketed(rng = z_rng, n_target = as.integer(n_contours))
  n_bins <- length(contour_breaks) - 1L

  bin_cols <- .make_diverging_bins_single_white(n_bins, z_rng, threshold, failure_dir)

  # Use the discrete contour bins for colors; map through numeric midpoints for steps legend
  bin_mid <- 0.5 * (contour_breaks[-1] + contour_breaks[-length(contour_breaks)])
  bin_vals <- scales::rescale(bin_mid, from = z_rng)

  # Ensure strictly increasing values for stepsn (avoid duplicates)
  eps <- 1e-9
  bin_vals <- cummax(bin_vals + seq_along(bin_vals) * eps)

  # Anchor ends so the scale covers the full legend range with no NA segment
  values_use <- c(0, pmin(pmax(bin_vals, 0), 1), 1)
  colors_use <- c(bin_cols[1], bin_cols, bin_cols[length(bin_cols)])

  # Remove any duplicate/unsorted values (defensive)
  ord <- order(values_use)
  values_use <- values_use[ord]
  colors_use <- colors_use[ord]

  # Collapse duplicates by keeping the LAST color for each value (prevents regularize.values warnings)
  keep <- !duplicated(values_use, fromLast = TRUE)
  values_use <- values_use[keep]
  colors_use <- colors_use[keep]

  # Final strict increasing nudge
  values_use <- cummax(values_use + seq_along(values_use) * eps)
  values_use <- pmin(values_use, 1)

  # ---------------------------------------------------------------------------
  # Legend breaks
  # ---------------------------------------------------------------------------
  # Rule of thumb: ~100 ticks or ~4 ticks per contour bin, whichever is larger, capped
  #n_ticks <- 7 #min(250L, max(25L, 4L * n_contours, 100L))
  #z_breaks_auto <- seq(z_rng[1], z_rng[2], length.out = n_ticks)

  legend_labeller <- .make_legend_labeller(max_labels = LEGEND_MAX_LABELS, white_value = threshold)

  # ---------------------------------------------------------------------------
  # Layout + theme
  # ---------------------------------------------------------------------------
  xy_ratio <- diff(range(x_breaks)) / diff(range(y_breaks))

  base_size <- 18 * text_size
  theme_surface <- theme_bw(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      plot.title = element_text(size = base_size + 3, hjust = 0),
      axis.title = element_text(size = base_size + 2),
      axis.text  = element_text(size = base_size - 1),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text  = element_text(size = base_size - 1),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(b = 6),
      plot.margin = margin(15, 25, 15, 15, unit = "pt"),
      legend.ticks = element_line(colour = "gray60", linewidth = 0.4),
      legend.ticks.length = grid::unit(2.5, "mm")
    )

  guide_fill <- guide_coloursteps(
    direction   = "horizontal",
    barwidth    = grid::unit(panel_size_in, "in"),
    barheight   = grid::unit(legend_barheight_in, "in"),
    show.limits = TRUE,
    ticks       = TRUE
  )

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_contour_filled(
      aes(z = .data[[z_var]], fill = after_stat(level_mid)),
      breaks = contour_breaks) +
    geom_contour(aes(z = .data[[z_var]]), breaks = threshold,
                 color = "black", linewidth = 1) +
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
      labels = legend_labeller,
      oob    = scales::squish,
      guide  = guide_fill
    ) +
    coord_fixed(ratio = xy_ratio, expand = FALSE) +
    labs(x = x_label, y = y_label, title = title, fill = "") +
    theme_surface

  k <- 1L
  if (isTRUE(facet)) {

    p <- p +
      facet_wrap(vars(!!sym(facet_by))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
    v <- data[[facet_by]]
    k <- if (is.factor(v)) nlevels(v) else length(unique(v))

  }

  dims <- .compute_facet_dims(
    k = k,
    panel_w_in = panel_in,
    panel_h_in = panel_in,
    left_right_in = left_right_in,
    top_in = top_in,
    bottom_in = bottom_in,
    gap_w_in = 0.25,
    gap_h_in = 0.25,
    max_cols = 2
  )

  list(p = p, width = dims$width, height = dims$height)
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
    size = 2.4,
    stroke = 1.2,
    show_legend = TRUE,
    ellipse_levels = NULL,
    ellipse_group = c("none", "scenario", "horizon", "scenario_horizon"),
    ellipse_color = "gray30",
    ellipse_linetype = "dashed",
    ellipse_linewidth = 0.5
) {

  ellipse_group <- match.arg(ellipse_group)


  panel_in = 6.5
  left_right_in <- 1.5
  bottom_in <- 1.0
  top_in <- 1.0

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

  # sanity: if factors present, enforce complete mappings
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
    labs(color = "GCM\nScenarios", shape = "Time\nHorizon")

  # ---------------------------------------------------------------------------
  # Optional ellipses
  # ---------------------------------------------------------------------------
  if (!is.null(ellipse_levels)) {
    if (!is.numeric(ellipse_levels) || any(!is.finite(ellipse_levels)) ||
        any(ellipse_levels <= 0 | ellipse_levels >= 1)) {
      stop("'ellipse_levels' must be numeric values strictly between 0 and 1.", call. = FALSE)
    }

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
      grp <- switch(
        ellipse_group,
        scenario = scenario_var,
        horizon = horizon_var,
        scenario_horizon = interaction(g[[scenario_var]], g[[horizon_var]], drop = TRUE)
      )

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
  }

  # ---------------------------------------------------------------------------
  # Legend control (keep surface legend from base; this only toggles GCM legends)
  # ---------------------------------------------------------------------------
  if (isTRUE(show_legend)) {
    p2 <- p2 + guides(
      color = guide_legend(order = 2, position = "right", direction = "vertical"),
      shape = guide_legend(order = 3, position = "right", direction = "vertical")
    )
  } else {
    p2 <- p2 + guides(color = "none", shape = "none")
  }

  list(p = p2, width = panel_in + left_right_in,
       height = panel_in + bottom_in + top_in)
}
