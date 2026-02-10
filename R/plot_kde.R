#' Plot a 2D KDE surface as a heatmap with contour lines
#'
#' @description
#' Plots a gridded 2D kernel density estimate (or any gridded surface) as a raster
#' heatmap with contour lines, faceted by a grouping variable (e.g., scenario).
#' Optionally overlays the original sample points.
#'
#' @details
#' Expected data shape:
#' - `kde_grid` is in long format: one row per grid point per facet.
#' - `samples` contains point locations to overlay (typically the ensemble points
#'   used to build the KDE).
#' - `x_col` and `y_col` must exist in both `kde_grid` and `samples`.
#' - `facet_col` and `z_col` must exist in `kde_grid`.
#'
#' Notes:
#' - Contours are computed from `z_col` within each facet.
#' - `fill_label` controls the legend title (via the fill scale).
#' - `x_limits` / `y_limits` use `coord_cartesian()` (zoom without dropping data).
#'
#' @param kde_grid Data frame in long format containing gridded surface values.
#' @param samples Data frame of sample points to overlay.
#' @param x_col Character scalar; x-axis column name in both inputs (e.g., temperature).
#' @param y_col Character scalar; y-axis column name in both inputs (e.g., precipitation).
#' @param facet_col Character scalar; column in `kde_grid` used for faceting (e.g., scenario).
#' @param z_col Character scalar; column in `kde_grid` used as the surface value (e.g., weight/density).
#' @param x_limits Numeric length-2 or `NULL`; x-axis limits passed to `coord_cartesian()`.
#' @param y_limits Numeric length-2 or `NULL`; y-axis limits passed to `coord_cartesian()`.
#' @param n_contours Integer >= 1; number of contour bins for `geom_contour()`.
#' @param interpolate_raster Logical; whether to interpolate raster pixels in `geom_raster()`.
#' @param points_alpha Numeric in [0, 1]; alpha for sample point overlay.
#' @param palette Character; viridis palette option passed to `scale_fill_viridis_c()`.
#' @param palette_begin Numeric in [0, 1]; start of viridis palette.
#' @param palette_end Numeric in [0, 1]; end of viridis palette.
#' @param fill_opacity Numeric in [0, 1]; alpha for the fill scale.
#' @param fill_label Character; legend title for the fill scale.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # Minimal reproducible example
#'
#' set.seed(1)
#'
#' # Sample points (e.g., ensemble projections)
#' samples <- data.frame(
#'   scenario = rep(c("SSP1", "SSP2"), each = 20),
#'   tavg = c(rnorm(20, 0, 1), rnorm(20, 1, 1)),
#'   prcp = c(rnorm(20, 0, 1), rnorm(20, 0.5, 1))
#' )
#'
#' # Simple evaluation grid
#' grid <- expand.grid(
#'   tavg = seq(-3, 3, length.out = 40),
#'   prcp = seq(-3, 3, length.out = 40),
#'   scenario = c("SSP1", "SSP2")
#' )
#'
#' # Fake KDE-like surface (for plotting only)
#' grid$weight <- with(
#'   grid,
#'   exp(-0.5 * (tavg^2 + prcp^2))
#' )
#'
#' # Plot
#' plot_kde(
#'   kde_grid = grid,
#'   samples = samples,
#'   x_col = "tavg",
#'   y_col = "prcp",
#'   facet_col = "scenario",
#'   z_col = "weight",
#'   n_contours = 10,
#'   interpolate_raster = FALSE,
#'   points_alpha = 0.6,
#'   fill_label = "KDE weight"
#' )
#' @import ggplot2 dplyr
#' @export
plot_kde <- function(
    kde_grid,
    samples,
    x_col = "tavg",
    y_col = "prcp",
    facet_col = "scenario",
    z_col = "weight",
    x_limits = NULL,
    y_limits = NULL,
    n_contours = 20,
    interpolate_raster = FALSE,
    points_alpha = 0.8,
    palette = "viridis",
    palette_begin = 0,
    palette_end = 1,
    fill_opacity = 0.8,
    fill_label = "Weight"
) {
  # ---------------------------------------------------------------------------
  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (!is.data.frame(kde_grid)) stop("kde_grid must be a data.frame.", call. = FALSE)
  if (!is.data.frame(samples)) stop("samples must be a data.frame.", call. = FALSE)

  req_kde <- c(x_col, y_col, facet_col, z_col)
  missing_kde <- setdiff(req_kde, names(kde_grid))
  if (length(missing_kde) > 0) {
    stop("kde_grid is missing required columns: ",
         paste(missing_kde, collapse = ", "), call. = FALSE)
  }

  req_samples <- c(x_col, y_col)
  missing_samples <- setdiff(req_samples, names(samples))
  if (length(missing_samples) > 0) {
    stop("samples is missing required columns: ",
         paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  .validate_limits <- function(lim, name) {
    if (is.null(lim)) return(NULL)
    if (!is.numeric(lim) || length(lim) != 2 || any(!is.finite(lim))) {
      stop(name, " must be NULL or numeric length-2 with finite values.", call. = FALSE)
    }
    if (lim[1] > lim[2]) stop(name, " must have lim[1] <= lim[2].", call. = FALSE)
    lim
  }

  x_limits <- .validate_limits(x_limits, "x_limits")
  y_limits <- .validate_limits(y_limits, "y_limits")

  if (!is.numeric(n_contours) || length(n_contours) != 1 || !is.finite(n_contours) || n_contours < 1) {
    stop("n_contours must be a single integer >= 1.", call. = FALSE)
  }
  n_contours <- as.integer(n_contours)

  if (!is.logical(interpolate_raster) || length(interpolate_raster) != 1) {
    stop("interpolate_raster must be a single logical.", call. = FALSE)
  }

  if (!is.numeric(points_alpha) || length(points_alpha) != 1 ||
      !is.finite(points_alpha) || points_alpha < 0 || points_alpha > 1) {
    stop("points_alpha must be a single numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(fill_opacity) || length(fill_opacity) != 1 ||
      !is.finite(fill_opacity) || fill_opacity < 0 || fill_opacity > 1) {
    stop("fill_opacity must be a single numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(palette_begin) || length(palette_begin) != 1 ||
      !is.finite(palette_begin) || palette_begin < 0 || palette_begin > 1) {
    stop("palette_begin must be a single numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(palette_end) || length(palette_end) != 1 ||
      !is.finite(palette_end) || palette_end < 0 || palette_end > 1) {
    stop("palette_end must be a single numeric in [0, 1].", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  ggplot(
    kde_grid,
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]]
    )
  ) +
    theme_bw() +
    geom_raster(
      aes(fill = .data[[z_col]]),
      interpolate = interpolate_raster
    ) +
    geom_contour(
      aes(z = .data[[z_col]]),
      bins = n_contours,
      linewidth = 0.3
    ) +
    geom_point(
      data = samples,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]]
      ),
      alpha = points_alpha,
      inherit.aes = FALSE
    ) +
    labs(
      x = "ΔT",
      y = "ΔP",
      fill = fill_label
    ) +
    facet_wrap(stats::as.formula(paste("~", facet_col))) +
    scale_fill_viridis_c(
      option = palette,
      begin = palette_begin,
      end = palette_end,
      alpha = fill_opacity,
      name = fill_label
    ) +
    coord_cartesian(
      expand = FALSE,
      xlim = x_limits,
      ylim = y_limits,
      clip = "on"
    )
}
