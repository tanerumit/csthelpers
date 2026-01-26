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
    legend_barwidth_in = 3.0,
    legend_barheight_in = 0.25,
    text_size = 0.7,
    facet = FALSE,
    facet_by = NULL,
    facet_levels = NULL
) {
  BASELINE_TOL <- 1e-10
  LEGEND_MAX_LABELS <- 14L

  COLOR_FAILURE <- "#df0000"
  COLOR_SAFE <- "#0033FF"
  COLOR_MID <- "#FFFFFF"
  COLOR_FAIL_LIGHT <- "#FEE5D9"
  COLOR_SAFE_LIGHT <- "#D6E3FF"

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

    nice_frac <- if (frac <= 1) 1 else if (frac <= 2) 2 else if (frac <= 2.5) 2.5 else if (frac <= 5) 5 else 10
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

  .make_diverging_bins_single_white <- function(n_bins, rng, thr, fail_dir) {
    lo <- rng[1]
    hi <- rng[2]
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo || n_bins <= 0) {
      return(rep(COLOR_MID, max(1L, n_bins)))
    }

    thr_c <- min(max(thr, lo), hi)
    prop <- (thr_c - lo) / (hi - lo)

    n_low <- max(1L, floor(n_bins * prop))
    n_high <- n_bins - n_low

    if (n_high < 1L) { n_high <- 1L; n_low <- n_bins - 1L }
    if (n_low < 1L)  { n_low  <- 1L; n_high <- n_bins - 1L }

    if (fail_dir == 1) {
      cols_low  <- if (n_low == 1L) COLOR_FAILURE else colorRampPalette(c(COLOR_FAILURE, COLOR_FAIL_LIGHT))(n_low)
      cols_high <- if (n_high == 1L) COLOR_MID else c(COLOR_MID, colorRampPalette(c(COLOR_SAFE_LIGHT, COLOR_SAFE))(n_high - 1L))
    } else {
      cols_low  <- if (n_low == 1L) COLOR_SAFE else colorRampPalette(c(COLOR_SAFE, COLOR_SAFE_LIGHT))(n_low)
      cols_high <- if (n_high == 1L) COLOR_MID else c(COLOR_MID, colorRampPalette(c(COLOR_FAIL_LIGHT, COLOR_FAILURE))(n_high - 1L))
    }

    cols <- c(cols_low, cols_high)
    if (length(cols) < n_bins) cols <- c(cols, rep(tail(cols, 1), n_bins - length(cols)))
    if (length(cols) > n_bins) cols <- cols[seq_len(n_bins)]
    cols
  }

  .make_legend_labeller <- function(max_labels = LEGEND_MAX_LABELS) {
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
      if (idx_keep[length(idx_keep)] != n) idx_keep <- c(idx_keep, n)

      pos <- which(finite)[idx_keep]
      out[pos] <- format(b[idx_keep], trim = TRUE, scientific = FALSE)
      out
    }
  }

  # ---------------------------------------------------------------------------
  # Validation
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

  z_rng <- if (is.null(z_limits)) z_rng_data else z_limits
  if (!is.null(z_limits)) data[[z_var]] <- pmin(pmax(data[[z_var]], z_limits[1]), z_limits[2])

  if (is.null(threshold)) {
    baseline_idx <- which(abs(data[[x_var]]) <= BASELINE_TOL & abs(data[[y_var]]) <= BASELINE_TOL)
    if (length(baseline_idx) == 0) {
      stop("Cannot infer threshold: no baseline point (x≈0,y≈0). Provide 'threshold'.", call. = FALSE)
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

  # ---------------------------------------------------------------------------
  # Theme + guide (will be OVERRIDDEN at save-time by ggsave_smart)
  # ---------------------------------------------------------------------------
  x_span <- diff(range(x_breaks))
  y_span <- diff(range(y_breaks))
  xy_ratio <- if (is.finite(x_span) && is.finite(y_span) && y_span > 0) x_span / y_span else 1

  base_size <- 12 * text_size

  theme_surface <- theme_bw(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      plot.title = element_text(size = base_size + 1.5, hjust = 0),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size - 1),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2.5),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(b = 6),
      plot.margin = margin(15, 25, 15, 15, unit = "pt"),
      legend.ticks = element_line(colour = "gray60", linewidth = 0.4),
      legend.ticks.length = grid::unit(2.5, "mm")
    )

  guide_fill <- guide_coloursteps(
    direction = "horizontal",
    barwidth  = grid::unit(legend_barwidth_in, "in"),
    barheight = grid::unit(legend_barheight_in, "in"),
    show.limits = TRUE,
    ticks = TRUE
  )

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

  if (isTRUE(facet)) {
    p <- p +
      facet_wrap(vars(!!sym(facet_by))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
  }

  # Metadata for save-time scaling
  attr(p, "legend_barwidth_in")  <- legend_barwidth_in
  attr(p, "legend_barheight_in") <- legend_barheight_in
  attr(p, "legend_max_labels") <- LEGEND_MAX_LABELS

  p
}


ggsave_smart <- function(
    filename,
    plot,
    target = c("journal", "ppt"),
    # Panel geometry (inches)
    panel_w_in = 2.2,
    panel_h_in = 2.0,
    gap_w_in = 0.18,
    gap_h_in = 0.18,
    margin_left_in = 0.7,
    margin_right_in = 0.7,
    margin_top_in = 0.7,
    margin_bottom_in = 0.7,
    # Journal sizing caps (inches)
    journal_width_in = c(single = 3.4, double = 6.7),
    journal_mode = c("double", "single"),
    # PPT slide size (inches)
    ppt_width_in = 13.33,
    ppt_height_in = 7.5,
    # Output
    dpi = NULL,
    bg = "white",
    device = NULL,
    verbose = FALSE
) {
  target <- match.arg(target)
  journal_mode <- match.arg(journal_mode)

  if (!is.character(filename) || length(filename) != 1 || nchar(filename) == 0) {
    stop("'filename' must be a non-empty character scalar.", call. = FALSE)
  }
  if (!inherits(plot, "ggplot")) stop("'plot' must be a ggplot object.", call. = FALSE)

  out_dir <- dirname(filename)
  if (!dir.exists(out_dir)) {
    ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) stop("Failed to create output directory: ", out_dir, call. = FALSE)
  }

  .is_num1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  .assert_num1_gt0 <- function(x, nm) {
    if (!.is_num1(x) || x <= 0) stop("'", nm, "' must be a single number > 0.", call. = FALSE)
  }

  .assert_num1_gt0(panel_w_in, "panel_w_in")
  .assert_num1_gt0(panel_h_in, "panel_h_in")
  .assert_num1_gt0(gap_w_in, "gap_w_in")
  .assert_num1_gt0(gap_h_in, "gap_h_in")
  .assert_num1_gt0(margin_left_in, "margin_left_in")
  .assert_num1_gt0(margin_right_in, "margin_right_in")
  .assert_num1_gt0(margin_top_in, "margin_top_in")
  .assert_num1_gt0(margin_bottom_in, "margin_bottom_in")

  if (!is.numeric(journal_width_in) || length(journal_width_in) != 2 || any(!is.finite(journal_width_in))) {
    stop("'journal_width_in' must be numeric length-2: c(single, double).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Infer facet layout / number of panels
  # ---------------------------------------------------------------------------
  .infer_layout <- function(p) {
    built <- ggplot2::ggplot_build(p)
    lay <- built$layout$layout
    n_panels <- if (!is.null(lay$PANEL)) length(unique(lay$PANEL)) else nrow(lay)
    n_panels <- max(1L, as.integer(n_panels))

    nrow <- 1L
    ncol <- 1L
    facet_params <- tryCatch(p$facet$params, error = function(e) NULL)

    if (!is.null(facet_params)) {
      if (!is.null(facet_params$nrow) && is.numeric(facet_params$nrow) && facet_params$nrow > 0) nrow <- as.integer(facet_params$nrow)
      if (!is.null(facet_params$ncol) && is.numeric(facet_params$ncol) && facet_params$ncol > 0) ncol <- as.integer(facet_params$ncol)
    }

    if (n_panels > 1L) {
      if (nrow > 1L && ncol == 1L) ncol <- as.integer(ceiling(n_panels / nrow))
      if (ncol > 1L && nrow == 1L) nrow <- as.integer(ceiling(n_panels / ncol))
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

  .compute_size <- function(nrow, ncol, panel_w, panel_h) {
    width <- margin_left_in + margin_right_in + ncol * panel_w + (ncol - 1L) * gap_w_in
    height <- margin_top_in + margin_bottom_in + nrow * panel_h + (nrow - 1L) * gap_h_in
    list(width = width, height = height)
  }

  size0 <- .compute_size(nrow, ncol, panel_w_in, panel_h_in)
  width_in <- size0$width
  height_in <- size0$height

  ext <- tolower(tools::file_ext(filename))

  if (is.null(device)) {
    device <- switch(
      ext,
      pdf  = grDevices::cairo_pdf,
      png  = "png",
      tiff = "tiff",
      tif  = "tiff",
      svg  = "svg",
      eps  = "eps",
      NULL
    )
  }

  if (target == "journal") {
    if (is.null(dpi)) dpi <- if (ext %in% c("png", "tif", "tiff")) 300 else NA_real_

    cap_w <- unname(journal_width_in[[journal_mode]])
    if (!.is_num1(cap_w) || cap_w <= 0) stop("Invalid journal width cap.", call. = FALSE)

    if (n_panels == 1L) {
      # Single panel: enforce publication width, scale height proportionally
      scale_fac <- cap_w / width_in
      width_in <- cap_w
      height_in <- height_in * scale_fac
      if (verbose) message("Single-panel journal: set width to ", cap_w, " in; scaled height.")
    } else if (width_in > cap_w) {
      fixed_w <- margin_left_in + margin_right_in + (ncol - 1L) * gap_w_in
      avail_w <- cap_w - fixed_w
      if (avail_w <= 0) stop("Margins/gaps exceed journal width cap; reduce margins or gaps.", call. = FALSE)

      scale_fac <- avail_w / (ncol * panel_w_in)
      panel_w_in2 <- panel_w_in * scale_fac
      panel_h_in2 <- panel_h_in * scale_fac

      size1 <- .compute_size(nrow, ncol, panel_w_in2, panel_h_in2)
      width_in <- size1$width
      height_in <- size1$height

      if (verbose) message("Faceted journal: scaled panels by factor ", signif(scale_fac, 3), ".")
    }
  } else {
    if (is.null(dpi)) dpi <- 300
    width_in <- ppt_width_in
    height_in <- ppt_height_in
  }

  if (!.is_num1(width_in) || !.is_num1(height_in) || width_in <= 0 || height_in <= 0) {
    stop("Computed invalid output size.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Panel-area dimensions (THIS is what you expect legend width to track)
  # ---------------------------------------------------------------------------
  panel_area_w_in <- width_in - margin_left_in - margin_right_in
  panel_area_h_in <- height_in - margin_top_in - margin_bottom_in

  # Defensive bounds
  panel_area_w_in <- max(1e-6, panel_area_w_in)
  panel_area_h_in <- max(1e-6, panel_area_h_in)

  # ---------------------------------------------------------------------------
  # Override legend bar size using plot attributes:
  # - if spec <= 1: fraction of PANEL AREA, not figure width
  # - if spec > 1: absolute inches
  # ---------------------------------------------------------------------------
  .spec_to_in <- function(spec, ref_in, default_in) {
    if (!is.numeric(spec) || length(spec) != 1 || !is.finite(spec) || spec <= 0) return(default_in)
    if (spec <= 1) return(spec * ref_in)
    spec
  }

  bw_spec <- attr(plot, "legend_barwidth_spec", exact = TRUE)
  bh_spec <- attr(plot, "legend_barheight_spec", exact = TRUE)

  barwidth_in <- .spec_to_in(bw_spec, panel_area_w_in, default_in = 0.85 * panel_area_w_in)
  barheight_in <- .spec_to_in(bh_spec, panel_area_h_in, default_in = 0.04 * panel_area_h_in)

  # Clamp: never exceed panel width; keep legible height
  barwidth_in <- max(1.5, min(panel_area_w_in, barwidth_in))
  barheight_in <- max(0.10, min(0.40, barheight_in))

  plot <- plot + ggplot2::guides(
    fill = ggplot2::guide_coloursteps(
      direction = "horizontal",
      barwidth  = grid::unit(barwidth_in, "in"),
      barheight = grid::unit(barheight_in, "in"),
      show.limits = TRUE,
      ticks = TRUE
    )
  )

  # ---------------------------------------------------------------------------
  # Save with cairo fallback for PDF
  # ---------------------------------------------------------------------------
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
    if (is.finite(dpi)) args$dpi <- dpi
    do.call(ggplot2::ggsave, args)
    TRUE
  }

  ok <- FALSE
  err <- NULL

  if (identical(ext, "pdf") && identical(device, grDevices::cairo_pdf)) {
    ok <- tryCatch(.ggsave_do(grDevices::cairo_pdf), error = function(e) { err <<- conditionMessage(e); FALSE })
    if (!ok) {
      if (verbose) message("cairo_pdf failed (", err, "); falling back to grDevices::pdf().")
      ok <- tryCatch(.ggsave_do(grDevices::pdf), error = function(e) { err <<- conditionMessage(e); FALSE })
    }
  } else {
    ok <- tryCatch(.ggsave_do(device), error = function(e) { err <<- conditionMessage(e); FALSE })
  }

  if (!ok) stop("ggsave failed: ", err, call. = FALSE)

  invisible(list(
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
