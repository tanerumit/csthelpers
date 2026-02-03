# =============================================================================
# Climate Response Surface Visualization Functions
# =============================================================================

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

#' Enforce strict monotonicity on a numeric vector
#' @keywords internal
.enforce_strict_monotonic <- function(x, eps = 1e-9) {
  n <- length(x)
  if (n < 2L) return(x)
  for (i in 2L:n) {
    if (x[i] <= x[i - 1L]) {
      x[i] <- x[i - 1L] + eps
    }
  }
  x
}


#' Generate "nice" bracketed breaks for contours
#' @keywords internal
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

  n_bins_target <- as.integer(max(1L, round(n_target)))
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

  dec <- max(0L, -floor(log10(step)) + 2L)
  lo <- round(floor(lo0 / step) * step, dec)
  hi <- round(ceiling(hi0 / step) * step, dec)

  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    stop("Internal error: failed to compute valid bracketed limits.", call. = FALSE)
  }

  br <- seq(lo, hi, by = step)
  br <- round(br, dec)

  br
}


#' Create diverging color bins centered on threshold
#' @keywords internal
.make_diverging_bins <- function(n_bins, rng, thr, fail_dir,
                                 col_failure = "#df0000",
                                 col_safe = "#0033FF",
                                 col_mid = "#FFFFFF",
                                 col_fail_light = "#FEE5D9",
                                 col_safe_light = "#D6E3FF") {
  lo <- rng[1]
  hi <- rng[2]

  if (!is.finite(lo) || !is.finite(hi) || hi <= lo || n_bins < 1L) {
    return(rep(col_mid, max(1L, n_bins)))
  }
  if (n_bins < 2L) {
    return(col_mid)
  }

  thr_c <- min(max(thr, lo), hi)
  prop <- (thr_c - lo) / (hi - lo)
  prop <- max(0.05, min(0.95, prop))

  n_low <- max(1L, floor(n_bins * prop))
  n_high <- n_bins - n_low

  # Ensure both have at least 1 bin
  if (n_high < 1L) {
    n_high <- 1L
    n_low <- max(1L, n_bins - 1L)
  }
  if (n_low < 1L) {
    n_low <- 1L
    n_high <- max(1L, n_bins - 1L)
  }

  if (fail_dir == 1) {
    # Safe below threshold, failure above
    cols_low <- if (n_low == 1L) {
      col_mid
    } else {
      ramp_safe <- colorRampPalette(c(col_safe, col_safe_light))
      c(ramp_safe(n_low - 1L), col_mid)
    }
    cols_high <- if (n_high == 1L) {
      col_failure
    } else {
      ramp_fail <- colorRampPalette(c(col_fail_light, col_failure))
      ramp_fail(n_high)
    }
  } else {
    # Failure below threshold, safe above
    cols_low <- if (n_low == 1L) {
      col_mid
    } else {
      ramp_fail <- colorRampPalette(c(col_failure, col_fail_light))
      c(ramp_fail(n_low - 1L), col_mid)
    }
    cols_high <- if (n_high == 1L) {
      col_safe
    } else {
      ramp_safe <- colorRampPalette(c(col_safe_light, col_safe))
      ramp_safe(n_high)
    }
  }

  cols <- c(cols_low, cols_high)

  # Ensure exact length
  if (length(cols) < n_bins) {
    cols <- c(cols, rep(cols[length(cols)], n_bins - length(cols)))
  } else if (length(cols) > n_bins) {
    cols <- cols[seq_len(n_bins)]
  }

  cols
}


#' Create legend labeller with maximum label count
#' @keywords internal
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

  # Final monotonicity enforcement and clamping
  values_use <- .enforce_strict_monotonic(values_use, eps = 1e-9)
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
  attr(p, "x_var")                 <- x_var
  attr(p, "y_var")                 <- y_var
  attr(p, "x_breaks")              <- x_breaks
  attr(p, "y_breaks")              <- y_breaks

  p
}


# -----------------------------------------------------------------------------
# GCM Overlay Helpers (REVISED)
# -----------------------------------------------------------------------------

#' Get group identifier for ellipse/KDE grouping
#' @keywords internal
.get_group_id <- function(df, mode, scenario_var, horizon_var) {
  if (mode == "none") return(rep("all", nrow(df)))
  if (mode == "scenario") return(as.character(df[[scenario_var]]))
  if (mode == "horizon") return(as.character(df[[horizon_var]]))
  # scenario_horizon
  paste0(as.character(df[[scenario_var]]), " / ", as.character(df[[horizon_var]]))
}

#' Compute ellipse coordinates from mean and covariance
#' @keywords internal
.ellipse_df_from_mu_S <- function(mu, S, level, n = 181L) {
  r <- sqrt(stats::qchisq(level, df = 2))
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  A <- ee$vectors %*% diag(sqrt(vals), 2, 2)
  t_seq <- seq(0, 2 * pi, length.out = n)
  circle <- rbind(cos(t_seq), sin(t_seq))
  pts <- t(mu + A %*% circle)
  data.frame(x = pts[, 1], y = pts[, 2], stringsAsFactors = FALSE)
}

#' Compute robust mean and covariance using MCD
#' @keywords internal
.robust_mu_S <- function(x, y) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for robust ellipses (MASS::cov.rob).", call. = FALSE)
  }
  m <- cbind(x, y)
  if (nrow(m) < 5L) return(NULL)

  fit <- tryCatch(
    MASS::cov.rob(m, method = "mcd"),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)
  if (any(!is.finite(fit$center)) || any(!is.finite(fit$cov))) return(NULL)

  # FIX 4.1: Check condition number for near-singular covariance
  cond_num <- tryCatch(kappa(fit$cov), error = function(e) Inf)
  if (!is.finite(cond_num) || cond_num > 1e10) return(NULL)

  list(mu = as.numeric(fit$center), S = fit$cov)
}

#' Compute standard mean and covariance
#' @keywords internal
.norm_mu_S <- function(x, y) {
  m <- cbind(x, y)
  if (nrow(m) < 3L) return(NULL)
  mu <- colMeans(m)
  S <- stats::cov(m)
  if (any(!is.finite(mu)) || any(!is.finite(S))) return(NULL)

  # FIX 4.1: Check condition number
  cond_num <- tryCatch(kappa(S), error = function(e) Inf)
  if (!is.finite(cond_num) || cond_num > 1e10) return(NULL)

  list(mu = as.numeric(mu), S = S)
}

#' Compute KDE HDR contours with plug-in or NRD bandwidth
#' @keywords internal
.kde_hdr_contours_df <- function(x, y,
                                 mass_levels,
                                 n = 120L,
                                 bw_method = c("auto", "plugin", "nrd"),
                                 bw_adjust = 1.0) {

  bw_method <- match.arg(bw_method)

  # Input validation and cleaning
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n_pts <- length(x)

  if (n_pts < 10L) return(NULL)

  # Standardize for numerical stability
  mx <- mean(x); my <- mean(y)
  sx <- sd(x); sy <- sd(y)

  if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) return(NULL)

  xs <- (x - mx) / sx
  ys <- (y - my) / sy

  data_std <- cbind(xs, ys)
  # Select bandwidth matrix on standardized scale
  H <- NULL
  if (bw_method == "plugin") {
    H_raw <- .kde_try_Hpi_2d(data_std)
    if (is.null(H_raw)) {
      warning("ks::Hpi failed; KDE contours not generated.", call. = FALSE)
      return(NULL)
    }
    H <- H_raw * (bw_adjust^2)
  } else if (bw_method == "nrd") {
    H_raw <- .kde_try_nrd_2d(xs, ys)
    if (is.null(H_raw)) {
      warning("MASS::bandwidth.nrd failed; KDE contours not generated.", call. = FALSE)
      return(NULL)
    }
    H <- H_raw * (bw_adjust^2)
  } else {
    # auto: try plugin first, fall back to nrd
    H <- .kde_pick_H_2d(data_std, bw_method = "auto", bw_adjust = bw_adjust)
    if (is.null(H)) return(NULL)
  }

  if (!.kde_is_pd_2x2(H)) return(NULL)

  # KDE estimation
  grid_x <- NULL
  grid_y <- NULL
  z <- NULL

  if (requireNamespace("ks", quietly = TRUE)) {
    tryCatch({
      kde_fit <- ks::kde(
        x = data_std,
        H = H,
        gridsize = c(n, n),
        xmin = c(min(xs) - 3 * sqrt(H[1, 1]), min(ys) - 3 * sqrt(H[2, 2])),
        xmax = c(max(xs) + 3 * sqrt(H[1, 1]), max(ys) + 3 * sqrt(H[2, 2]))
      )

      grid_x <- kde_fit$eval.points[[1]]
      grid_y <- kde_fit$eval.points[[2]]
      z <- kde_fit$estimate
    }, error = function(e) NULL)
  }

  # Fallback to MASS::kde2d if ks failed
  if (is.null(z)) {
    if (!requireNamespace("MASS", quietly = TRUE)) return(NULL)

    tryCatch({
      h_vec <- sqrt(diag(H))

      kd_mass <- MASS::kde2d(
        x = xs,
        y = ys,
        n = n,
        h = h_vec,
        lims = c(
          min(xs) - 3 * h_vec[1], max(xs) + 3 * h_vec[1],
          min(ys) - 3 * h_vec[2], max(ys) + 3 * h_vec[2]
        )
      )

      grid_x <- kd_mass$x
      grid_y <- kd_mass$y
      z <- kd_mass$z
    }, error = function(e) NULL)
  }

  if (is.null(z) || !all(is.finite(z))) return(NULL)

  # HDR threshold computation
  dx <- diff(grid_x[1:2])
  dy <- diff(grid_y[1:2])

  if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) return(NULL)

  z_vec <- as.vector(z)
  total_mass <- sum(z_vec) * dx * dy

  if (!is.finite(total_mass) || total_mass <= 0) return(NULL)

  # Normalize and compute cumulative mass
  z_norm <- z_vec / total_mass
  ord <- order(z_norm, decreasing = TRUE)
  z_sorted <- z_norm[ord]
  cum_mass <- cumsum(z_sorted * dx * dy)

  # Find density threshold for each mass level
  thresholds <- vapply(mass_levels, function(m) {
    idx <- which(cum_mass >= m)[1]
    if (is.na(idx)) min(z_sorted) else z_sorted[idx]
  }, numeric(1))

  # Convert back to unnormalized scale
  thresholds_raw <- thresholds * total_mass

  # FIX 3.4: Pre-allocate contour list
  cl_all <- vector("list", length(mass_levels) * 10L)
  idx <- 0L

  for (i in seq_along(thresholds_raw)) {
    lev <- thresholds_raw[i]

    # FIX 4.2: Check for NULL from contourLines
    cl <- tryCatch(
      grDevices::contourLines(grid_x, grid_y, z, levels = lev),
      error = function(e) NULL
    )

    if (is.null(cl) || length(cl) == 0L) next

    for (j in seq_along(cl)) {
      x_orig <- cl[[j]]$x * sx + mx
      y_orig <- cl[[j]]$y * sy + my

      idx <- idx + 1L
      cl_all[[idx]] <- data.frame(
        x = x_orig,
        y = y_orig,
        level = mass_levels[i],
        piece = paste0("hdr_", i, "_", j),
        stringsAsFactors = FALSE
      )
    }
  }

  if (idx == 0L) return(NULL)

  dplyr::bind_rows(cl_all[seq_len(idx)])
}


#' Overlay GCM information on a climate response surface plot
#'
#' @description
#' Adds Global Climate Model (GCM) summary information---typically scenario--horizon
#' change signals---to an existing climate response surface created with
#' \code{\link{climate_surface_base}}.
#'
#' The function overlays:
#' \itemize{
#'   \item GCM points positioned in the same \code{(x, y)} climate-change space
#'         as the response surface (e.g. \eqn{\Delta}precipitation vs
#'         \eqn{\Delta}temperature),
#'   \item scenario-specific colors and horizon-specific point shapes,
#'   \item optional spread summaries via ellipses or KDE highest density regions.
#' }
#'
#' The underlying response surface is \strong{not modified}; this function only adds
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
#' \code{gcm_data} before plotting, ensuring GCM points and spread summaries
#' are confined to the visible plotting domain.
#'
#' \strong{Spread methods:}
#' \itemize{
#'   \item \code{"none"}: No spread visualization.
#'   \item \code{"ellipse_norm"}: Standard Gaussian ellipses via \code{stat_ellipse()}.
#'   \item \code{"ellipse_robust"}: Robust ellipses using Minimum Covariance
#'         Determinant (MCD) via \code{MASS::cov.rob()}, less sensitive to outliers.
#'   \item \code{"kde"}: Kernel density estimation highest density regions (HDR)
#'         with optional plug-in bandwidth selection via \code{ks::Hpi()}.
#' }
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
#' Transparency of GCM points (0--1). Useful for dense ensembles.
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
#' @param spread_method
#' Method for visualizing ensemble spread:
#' \itemize{
#'   \item \code{"none"}: No spread visualization (default).
#'   \item \code{"ellipse_norm"}: Standard Gaussian ellipses.
#'   \item \code{"ellipse_robust"}: Robust (MCD-based) ellipses.
#'   \item \code{"kde"}: KDE highest density region contours.
#' }
#'
#' @param spread_levels
#' Numeric vector of probability levels (values strictly between 0 and 1)
#' for spread visualization. For ellipses, these are confidence levels;
#' for KDE, these are probability mass levels. Example: \code{c(0.5, 0.9)}.
#'
#' @param spread_group
#' Controls grouping for spread estimation:
#' \itemize{
#'   \item \code{"none"}: One spread summary using all GCM points.
#'   \item \code{"scenario"}: Separate summaries per scenario.
#'   \item \code{"horizon"}: Separate summaries per horizon.
#'   \item \code{"scenario_horizon"}: Separate summaries for each
#'         scenario--horizon combination.
#' }
#'
#' @param spread_color
#' Line color for spread visualization.
#'
#' @param spread_linetype
#' Line type for spread visualization (e.g. \code{"dashed"}, \code{"solid"}).
#'
#' @param spread_linewidth
#' Line width for spread visualization.
#'
#' @param kde_bw_method
#' Bandwidth selection method for KDE (only used when \code{spread_method = "kde"}):
#' \itemize{
#'   \item \code{"auto"}: Try plug-in first, fall back to NRD.
#'   \item \code{"plugin"}: Full 2D plug-in bandwidth via \code{ks::Hpi()}.
#'   \item \code{"nrd"}: Normal reference distribution via \code{MASS::bandwidth.nrd()}.
#' }
#'
#' @param kde_n
#' Grid resolution for KDE (integer >= 40). Only used when
#' \code{spread_method = "kde"}.
#'
#' @param kde_bw_adjust
#' Bandwidth adjustment factor for KDE (multiplier > 0). Only used when
#' \code{spread_method = "kde"}.
#'
#' @return
#' A ggplot object consisting of the original climate response surface with
#' GCM points and optional spread summaries overlaid.
#'
#' @seealso
#' \code{\link{climate_surface_base}}, \code{\link[ggplot2]{stat_ellipse}}
#'
#' @examples
#' \dontrun{
#' # Basic overlay with points only
#' climate_surface_gcm_overlay(p, gcm_data)
#'
#' # With Gaussian ellipses at 50% and 90% confidence
#' climate_surface_gcm_overlay(p, gcm_data,
#'   spread_method = "ellipse_norm",
#'   spread_levels = c(0.5, 0.9))
#'
#' # With robust ellipses grouped by scenario
#' climate_surface_gcm_overlay(p, gcm_data,
#'   spread_method = "ellipse_robust",
#'   spread_levels = 0.9,
#'   spread_group = "scenario")
#'
#' # With KDE HDR contours
#' climate_surface_gcm_overlay(p, gcm_data,
#'   spread_method = "kde",
#'   spread_levels = c(0.5, 0.9),
#'   kde_bw_method = "plugin")
#' }
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
    # Unified spread parameters
    spread_method = c("none", "ellipse_norm", "ellipse_robust", "kde"),
    spread_levels = NULL,
    spread_group = c("none", "scenario", "horizon", "scenario_horizon"),
    spread_color = "gray30",
    spread_linetype = "dashed",
    spread_linewidth = 0.5,
    # KDE-specific parameters
    kde_bw_method = c("auto", "plugin", "nrd"),
    kde_n = 120L,
    kde_bw_adjust = 1.0
) {

  # Match arguments
  spread_method <- match.arg(spread_method)
  spread_group <- match.arg(spread_group)
  kde_bw_method <- match.arg(kde_bw_method)

  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------

  # Check plot object first
  if (is.null(p) || !inherits(p, "ggplot")) {
    stop("'p' must be a ggplot object.", call. = FALSE)
  }

  # Early return for NULL data
  if (is.null(gcm_data)) return(p)

  if (!is.data.frame(gcm_data)) {
    stop("'gcm_data' must be a data.frame.", call. = FALSE)
  }

  # Infer x_var/y_var from plot attributes if not provided
  if (is.null(x_var)) x_var <- attr(p, "x_var")
  if (is.null(y_var)) y_var <- attr(p, "y_var")

  if (is.null(x_var) || is.null(y_var)) {
    stop("Provide 'x_var' and 'y_var' (or pass a plot from climate_surface_base() with these attributes).", call. = FALSE)
  }

  # Check required columns exist
  req <- c(x_var, y_var, scenario_var, horizon_var)
  miss <- setdiff(req, names(gcm_data))
  if (length(miss) > 0L) {
    stop("Missing columns in gcm_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # Check x/y are numeric
  if (!is.numeric(gcm_data[[x_var]]) || !is.numeric(gcm_data[[y_var]])) {
    stop("gcm_data x/y columns must be numeric.", call. = FALSE)
  }

  # Validate scalar parameters
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single number in [0,1].", call. = FALSE)
  }

  # Validate spread_levels if spread visualization is requested
  if (spread_method != "none" && !is.null(spread_levels)) {
    if (!is.numeric(spread_levels) || any(!is.finite(spread_levels)) ||
        any(spread_levels <= 0 | spread_levels >= 1)) {
      stop("'spread_levels' must be numeric values strictly between 0 and 1.", call. = FALSE)
    }
  }

  # Warn if spread_method is set but no levels provided
  if (spread_method != "none" && is.null(spread_levels)) {
    warning("'spread_method' is set but 'spread_levels' is NULL; no spread visualization will be drawn.",
            call. = FALSE)
  }

  # Validate KDE-specific parameters
  if (!is.numeric(kde_n) || length(kde_n) != 1L || !is.finite(kde_n) || kde_n < 40L) {
    stop("'kde_n' must be a single integer >= 40.", call. = FALSE)
  }
  kde_n <- as.integer(kde_n)

  if (!is.numeric(kde_bw_adjust) || length(kde_bw_adjust) != 1L ||
      !is.finite(kde_bw_adjust) || kde_bw_adjust <= 0) {
    stop("'kde_bw_adjust' must be a single number > 0.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Data preparation
  # ---------------------------------------------------------------------------

  # Filter for finite coordinates first
  g <- gcm_data %>%
    dplyr::filter(is.finite(.data[[x_var]]), is.finite(.data[[y_var]]))

  # Early return if no valid data
  if (nrow(g) == 0L) return(p)

  # Clip to base plot range when available
  x_breaks <- attr(p, "x_breaks")
  y_breaks <- attr(p, "y_breaks")

  if (is.numeric(x_breaks) && length(x_breaks) >= 2L) {
    g <- g %>% dplyr::filter(
      .data[[x_var]] >= min(x_breaks),
      .data[[x_var]] <= max(x_breaks)
    )
  }
  if (is.numeric(y_breaks) && length(y_breaks) >= 2L) {
    g <- g %>% dplyr::filter(
      .data[[y_var]] >= min(y_breaks),
      .data[[y_var]] <= max(y_breaks)
    )
  }

  if (nrow(g) == 0L) return(p)

  # Apply scenario level filtering and factorization
  if (!is.null(scenario_levels)) {
    g <- g %>% dplyr::filter(.data[[scenario_var]] %in% scenario_levels)
    if (nrow(g) == 0L) return(p)
    g[[scenario_var]] <- factor(g[[scenario_var]], levels = scenario_levels)
  }

  # Apply horizon level filtering and factorization
  if (!is.null(horizon_levels)) {
    g <- g %>% dplyr::filter(.data[[horizon_var]] %in% horizon_levels)
    if (nrow(g) == 0L) return(p)
    g[[horizon_var]] <- factor(g[[horizon_var]], levels = horizon_levels)
  }

  # Validate color/shape mappings AFTER factor creation
  if (is.factor(g[[scenario_var]])) {
    need <- levels(g[[scenario_var]])
    miss_cols <- setdiff(need, names(scenario_colors))
    if (length(miss_cols) > 0L) {
      stop("scenario_colors missing: ", paste(miss_cols, collapse = ", "), call. = FALSE)
    }
  }
  if (is.factor(g[[horizon_var]])) {
    need <- levels(g[[horizon_var]])
    miss_shapes <- setdiff(need, names(horizon_shapes))
    if (length(miss_shapes) > 0L) {
      stop("horizon_shapes missing: ", paste(miss_shapes, collapse = ", "), call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Add GCM points
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
  # Spread visualization (ellipses or KDE)
  # ---------------------------------------------------------------------------
  if (spread_method != "none" && !is.null(spread_levels)) {

    # Prepare grouping
    grp <- .get_group_id(g, spread_group, scenario_var, horizon_var)
    g2 <- g
    g2[[".spread_group"]] <- grp

    if (spread_method == "ellipse_norm") {
      # Standard Gaussian ellipses via stat_ellipse
      if (spread_group == "none") {
        for (lv in spread_levels) {
          p2 <- p2 + stat_ellipse(
            data = g,
            mapping = aes(x = .data[[x_var]], y = .data[[y_var]]),
            level = lv,
            type = "norm",
            color = spread_color,
            linetype = spread_linetype,
            linewidth = spread_linewidth
          )
        }
      } else {
        for (lv in spread_levels) {
          p2 <- p2 + stat_ellipse(
            data = g2,
            mapping = aes(
              x = .data[[x_var]],
              y = .data[[y_var]],
              group = .data[[".spread_group"]]
            ),
            level = lv,
            type = "norm",
            color = spread_color,
            linetype = spread_linetype,
            linewidth = spread_linewidth
          )
        }
      }

    } else if (spread_method == "ellipse_robust") {
      # Robust ellipses via MCD
      g2_split <- split(g2, g2[[".spread_group"]])
      n_groups <- length(g2_split)
      n_levels <- length(spread_levels)

      ell_list <- vector("list", n_groups * n_levels)
      idx <- 0L

      for (grp_name in names(g2_split)) {
        sub <- g2_split[[grp_name]]
        fit <- .robust_mu_S(sub[[x_var]], sub[[y_var]])
        if (is.null(fit)) next

        for (lv in spread_levels) {
          df <- .ellipse_df_from_mu_S(fit$mu, fit$S, level = lv, n = 181L)
          names(df) <- c("x", "y")
          df[[".spread_group"]] <- grp_name
          df[["level"]] <- lv
          df[["piece"]] <- paste0("ell_", grp_name, "_", format(lv))

          idx <- idx + 1L
          ell_list[[idx]] <- df
        }
      }

      if (idx > 0L) {
        ell_df <- dplyr::bind_rows(ell_list[seq_len(idx)])
        p2 <- p2 + geom_path(
          data = ell_df,
          mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["piece"]]),
          inherit.aes = FALSE,
          color = spread_color,
          linetype = spread_linetype,
          linewidth = spread_linewidth
        )
      }

    } else if (spread_method == "kde") {
      # KDE HDR contours
      g2_split <- split(g2, g2[[".spread_group"]])

      kde_list <- vector("list", length(g2_split))
      idx <- 0L

      for (grp_name in names(g2_split)) {
        sub <- g2_split[[grp_name]]

        d <- .kde_hdr_contours_df(
          x = sub[[x_var]],
          y = sub[[y_var]],
          mass_levels = sort(unique(spread_levels)),
          n = kde_n,
          bw_method = kde_bw_method,
          bw_adjust = kde_bw_adjust
        )

        if (is.null(d) || nrow(d) == 0L) next

        d[[".spread_group"]] <- grp_name
        idx <- idx + 1L
        kde_list[[idx]] <- d
      }

      if (idx > 0L) {
        kde_df <- dplyr::bind_rows(kde_list[seq_len(idx)])

        p2 <- p2 + geom_path(
          data = kde_df,
          mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["piece"]]),
          inherit.aes = FALSE,
          color = spread_color,
          linetype = spread_linetype,
          linewidth = spread_linewidth
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Legend control
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

