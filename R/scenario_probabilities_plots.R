

plot_kde_weights <- function(
    kde_long,
    obs_data,
    ta_col = "tavg",
    pr_col = "prcp",
    scenario_col = "scenario",
    value_col = "weight",
    xlim = NULL,
    ylim = NULL,
    bins = 20,
    raster_interpolate = FALSE,
    point_alpha = 0.8,
    fill_option = "viridis",
    fill_begin = 0,
    fill_end = 1,
    fill_alpha = 0.8,
    fill_name = "Weight"
) {
  # ---------------------------------------------------------------------------
  # Input validation (fail early)
  # ---------------------------------------------------------------------------
  if (!is.data.frame(kde_long)) stop("kde_long must be a data.frame.", call. = FALSE)

  req_kde <- c(ta_col, pr_col, scenario_col, value_col)
  missing_kde <- setdiff(req_kde, names(kde_long))
  if (length(missing_kde) > 0) {
    stop("kde_long is missing required columns: ",
         paste(missing_kde, collapse = ", "), call. = FALSE)
  }

  req_obs <- c(ta_col, pr_col)
  missing_obs <- setdiff(req_obs, names(obs_data))
  if (length(missing_obs) > 0) {
    stop("obs_data is missing required columns: ",
         paste(missing_obs, collapse = ", "), call. = FALSE)
  }

  .validate_lim <- function(lim, name) {
    if (is.null(lim)) return(NULL)
    if (!is.numeric(lim) || length(lim) != 2 || any(!is.finite(lim))) {
      stop(name, " must be NULL or numeric length-2 with finite values.", call. = FALSE)
    }
    if (lim[1] > lim[2]) stop(name, " must have lim[1] <= lim[2].", call. = FALSE)
    lim
  }

  xlim <- .validate_lim(xlim, "xlim")
  ylim <- .validate_lim(ylim, "ylim")

  if (!is.numeric(bins) || length(bins) != 1 || !is.finite(bins) || bins < 1) {
    stop("bins must be a single integer >= 1.", call. = FALSE)
  }
  bins <- as.integer(bins)

  if (!is.logical(raster_interpolate) || length(raster_interpolate) != 1) {
    stop("raster_interpolate must be a single logical.", call. = FALSE)
  }

  if (!is.numeric(point_alpha) || length(point_alpha) != 1 ||
      !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be a single numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(fill_alpha) || length(fill_alpha) != 1 ||
      !is.finite(fill_alpha) || fill_alpha < 0 || fill_alpha > 1) {
    stop("fill_alpha must be a single numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(fill_begin) || length(fill_begin) != 1 ||
      !is.finite(fill_begin) || fill_begin < 0 || fill_begin > 1) {
    stop("fill_begin must be a single numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(fill_end) || length(fill_end) != 1 ||
      !is.finite(fill_end) || fill_end < 0 || fill_end > 1) {
    stop("fill_end must be a single numeric in [0, 1].", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  ggplot2::ggplot(
    kde_long,
    ggplot2::aes(
      x = .data[[ta_col]],
      y = .data[[pr_col]]
    )
  ) +
    ggplot2::theme_light() +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[[value_col]]),
      interpolate = raster_interpolate
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data[[value_col]]),
      bins = bins,
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = obs_data,
      ggplot2::aes(
        x = .data[[ta_col]],
        y = .data[[pr_col]]
      ),
      alpha = point_alpha,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = "ΔT",
      y = "ΔP",
      fill = "weight"
    ) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", scenario_col))) +
    ggplot2::scale_fill_viridis_c(
      option = fill_option,
      begin = fill_begin,
      end = fill_end,
      alpha = fill_alpha,
      name = fill_name
    ) +
    ggplot2::coord_cartesian(
      expand = FALSE,
      xlim = xlim,
      ylim = ylim,
      clip = "on"
    )
}

