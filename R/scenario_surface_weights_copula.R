# ==============================================================================
# Climate Scenario Surface Weight Estimation — Copula (Gaussian; KDE marginals)
# ==============================================================================

# ------------------------------------------------------------------------------
# Copula-specific helpers
# ------------------------------------------------------------------------------

#' Validate a 2D bandwidth vector (marginal KDE bandwidths)
#' @keywords internal
.scenwgt_validate_bw_vec_1d <- function(bw_vec) {
  is.numeric(bw_vec) && length(bw_vec) == 2L &&
    all(is.finite(bw_vec)) && all(bw_vec > 0)
}

#' Pick 1D KDE bandwidth (base R) with guards and clamping
#'
#' Bandwidth is in the same units as x. Uses deterministic rules only.
#'
#' @keywords internal
.scenwgt_pick_bw_1d <- function(x, bw_method = c("nrd0", "nrd", "sj"),
                               bw_min = 0, bw_max = Inf) {

  bw_method <- match.arg(bw_method)
  x <- x[is.finite(x)]

  if (length(x) < 2L) return(NA_real_)

  bw_raw <- tryCatch({
    if (bw_method == "nrd0") {
      stats::bw.nrd0(x)
    } else if (bw_method == "nrd") {
      stats::bw.nrd(x)
    } else {
      stats::bw.SJ(x)
    }
  }, error = function(e) NA_real_)

  if (!is.finite(bw_raw) || bw_raw <= .Machine$double.eps) return(NA_real_)

  bw_use <- max(bw_raw, bw_min)
  bw_use <- min(bw_use, bw_max)

  if (!is.finite(bw_use) || bw_use <= .Machine$double.eps) return(NA_real_)
  bw_use
}

#' Weighted 1D KDE PDF (Gaussian kernel) evaluated at x_eval
#'
#' Complexity: O(length(x_eval) * length(x_obs)) time, O(chunk * n_obs) memory.
#'
#' @keywords internal
.scenwgt_kde1_pdf <- function(x_eval, x_obs, w_obs, bw, chunk_size = 5000L) {

  n_eval <- length(x_eval)
  out <- rep(NA_real_, n_eval)

  inv_bw <- 1 / bw
  norm_const <- inv_bw

  for (s in seq.int(1L, n_eval, by = chunk_size)) {
    e <- as.integer(min(s + chunk_size - 1L, n_eval))
    x_chunk <- x_eval[s:e]
    u <- outer(x_chunk, x_obs, `-`) * inv_bw
    k <- stats::dnorm(u)
    out[s:e] <- as.numeric(k %*% w_obs) * norm_const
  }

  out
}

#' Weighted 1D KDE CDF (Gaussian kernel) evaluated at x_eval
#'
#' Complexity: O(length(x_eval) * length(x_obs)) time, O(chunk * n_obs) memory.
#'
#' @keywords internal
.scenwgt_kde1_cdf <- function(x_eval, x_obs, w_obs, bw, chunk_size = 5000L) {

  n_eval <- length(x_eval)
  out <- rep(NA_real_, n_eval)

  inv_bw <- 1 / bw

  for (s in seq.int(1L, n_eval, by = chunk_size)) {
    e <- as.integer(min(s + chunk_size - 1L, n_eval))
    x_chunk <- x_eval[s:e]
    u <- outer(x_chunk, x_obs, `-`) * inv_bw
    k <- stats::pnorm(u)
    out[s:e] <- as.numeric(k %*% w_obs)
  }

  out
}

#' Weighted correlation (Pearson) with guards
#' @keywords internal
.scenwgt_wcor <- function(x, y, w) {

  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(x) < 3L) return(NA_real_)

  s <- sum(w)
  if (!is.finite(s) || s <= .Machine$double.eps) return(NA_real_)
  w <- w / s

  mx <- sum(w * x)
  my <- sum(w * y)

  dx <- x - mx
  dy <- y - my

  vx <- sum(w * dx * dx)
  vy <- sum(w * dy * dy)
  if (!is.finite(vx) || !is.finite(vy) || vx <= .Machine$double.eps || vy <= .Machine$double.eps) {
    return(NA_real_)
  }

  cxy <- sum(w * dx * dy)

  rho <- cxy / sqrt(vx * vy)
  if (!is.finite(rho)) return(NA_real_)

  # Guard against exactly +/-1 (singularity)
  rho <- max(min(rho, 0.9999), -0.9999)
  rho
}

#' Gaussian copula density (2D) from Gaussian scores z1,z2 and correlation rho
#'
#' Uses a stable log-density formulation.
#'
#' @keywords internal
.scenwgt_gaussian_copula_density <- function(z1, z2, rho) {

  if (!is.finite(rho)) return(rep(NA_real_, length(z1)))

  one_minus_r2 <- 1 - rho * rho
  if (!is.finite(one_minus_r2) || one_minus_r2 <= .Machine$double.eps) {
    return(rep(NA_real_, length(z1)))
  }

  q <- z1 * z1 - 2 * rho * z1 * z2 + z2 * z2
  log_c <- -0.5 * log(one_minus_r2) - 0.5 * (q / one_minus_r2) + 0.5 * (z1 * z1 + z2 * z2)

  exp(log_c)
}

#' Clamp pseudo-observations away from {0,1} for qnorm stability
#' @keywords internal
.scenwgt_clamp_unit <- function(u, eps = 1e-12) {
  u <- pmax(u, eps)
  u <- pmin(u, 1 - eps)
  u
}

# ------------------------------------------------------------------------------
# Main Copula function
# ------------------------------------------------------------------------------

#' Calculate Climate-Informed Scenario Surface Weights (Gaussian copula; KDE marginals)
#'
#' @description
#' For each group (e.g., SSP), fits a Gaussian copula to the dependence structure
#' between temperature and precipitation after transforming each marginal using a
#' weighted 1D Gaussian-kernel KDE CDF. The joint density is evaluated on
#' `scenario_grid` (evaluation surface). If `normalize=TRUE`, the resulting surface
#' weights are normalized to sum to 1 (optionally area-weighted for regular grids).
#'
#' Interpretation: outputs are *relative surface weights* induced by the fitted model
#' and the evaluation surface; they are not calibrated probabilities.
#'
#' @param ensemble_data Data frame containing ensemble projections.
#' @param scenario_grid Data frame defining the evaluation surface (not necessarily regular).
#' @param pr_col Column name for precipitation variable.
#' @param ta_col Column name for temperature variable.
#' @param group_col Column name for grouping (e.g., scenario).
#' @param bw Optional marginal bandwidth vector c(bw_ta, bw_pr) in original units.
#' @param bw_method 1D bandwidth method: "nrd0", "nrd", or "sj".
#' @param k Grid-step multipliers for bandwidth floor (same semantics as KDE).
#' @param alpha Multiplier for data-driven (1D) bandwidth.
#' @param bw_min Minimum bandwidth bounds (per marginal).
#' @param bw_max Maximum bandwidth bounds (per marginal).
#' @param min_samples Minimum observations required per group.
#' @param normalize Logical; normalize weights to sum to 1?
#' @param area_weight Area weighting method: "regular" or "none".
#' @param scale Scaling method: "none", "global", or "by_group".
#' @param weights_col Optional column with per-row weights.
#' @param support Optional list with ta and/or pr bounds.
#' @param chunk_size Chunk size used in marginal evaluations.
#' @param diagnostics Logical; attach diagnostic attributes?
#' @param verbose Logical; print progress messages?
#'
#' @return Data frame with scenario_grid columns plus one weight column per group.
#'
#' @export
compute_scenario_surface_weights_copula <- function(
    ensemble_data,
    scenario_grid,
    pr_col = "prcp",
    ta_col = "tavg",
    group_col = "scenario",
    bw = NULL,
    bw_method = c("nrd0", "nrd", "sj"),
    k = c(1.5, 2.0),
    alpha = 1.0,
    bw_min = c(0, 0),
    bw_max = c(Inf, Inf),
    min_samples = 5L,
    normalize = TRUE,
    area_weight = c("regular", "none"),
    scale = c("none", "global", "by_group"),
    weights_col = NULL,
    support = NULL,
    chunk_size = 5000L,
    diagnostics = FALSE,
    verbose = TRUE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  required_cols <- c(ta_col, pr_col, group_col)
  missing_cols <- setdiff(required_cols, names(ensemble_data))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in ensemble_data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!all(c(ta_col, pr_col) %in% names(scenario_grid))) {
    stop("Both temperature and precipitation columns must exist in scenario_grid.", call. = FALSE)
  }

  if (nrow(scenario_grid) == 0L) stop("scenario_grid cannot be empty.", call. = FALSE)

  bw_method <- match.arg(bw_method)
  area_weight <- match.arg(area_weight)
  scale <- match.arg(scale)

  if (!is.logical(normalize) || length(normalize) != 1L) stop("normalize must be a single logical.", call. = FALSE)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L) stop("diagnostics must be a single logical.", call. = FALSE)

  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 3L) {
    stop("min_samples must be a single integer >= 3.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  if (!is.null(bw) && !.scenwgt_validate_bw_vec_1d(bw)) {
    stop("If provided, bw must be numeric length-2 with finite values > 0: c(bw_ta, bw_pr).", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 2L || any(!is.finite(k)) || any(k <= 0)) {
    stop("k must be numeric length-2 with finite values > 0.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a single finite numeric value > 0.", call. = FALSE)
  }

  if (!is.numeric(bw_min) || length(bw_min) != 2L || any(!is.finite(bw_min)) || any(bw_min < 0)) {
    stop("bw_min must be numeric length-2 with finite values >= 0.", call. = FALSE)
  }

  if (!is.numeric(bw_max) || length(bw_max) != 2L || any(is.na(bw_max)) || any(is.nan(bw_max))) {
    stop("bw_max must be numeric length-2 with non-missing values.", call. = FALSE)
  }
  if (any(bw_max <= 0)) stop("bw_max must be > 0 (use Inf for no upper bound).", call. = FALSE)
  if (any(bw_max <= bw_min)) stop("bw_max must be greater than bw_min for both dimensions.", call. = FALSE)

  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || chunk_size < 100L) {
    stop("chunk_size must be a single integer >= 100.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  if (!is.null(weights_col)) {
    if (!is.character(weights_col) || length(weights_col) != 1L || !nzchar(weights_col)) {
      stop("weights_col must be NULL or a single non-empty character column name.", call. = FALSE)
    }
    if (!(weights_col %in% names(ensemble_data))) {
      stop("weights_col not found in ensemble_data: ", weights_col, call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Grid context (shared)
  # ---------------------------------------------------------------------------

  grid_ctx <- .scenwgt_prepare_grid_context(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    ta_col = ta_col,
    pr_col = pr_col,
    scale = scale,
    area_weight = area_weight,
    support = support,
    verbose = verbose,
    tag = "COPULA"
  )

  n_grid <- grid_ctx$n_grid
  support_mask <- grid_ctx$support_mask

  # ---------------------------------------------------------------------------
  # Loop over groups (split() for efficiency)
  # ---------------------------------------------------------------------------

  out <- scenario_grid
  existing_cols <- names(out)
  skipped <- character(0)

  diag_bw <- list()
  diag_rho <- list()
  diag_eff_n <- list()
  diag_scaling <- list()
  diag_bw_method <- list()

  data_split <- split(ensemble_data, ensemble_data[[group_col]])
  groups <- sort(names(data_split))

  for (grp in groups) {

    df <- data_split[[grp]]

    ta_full <- df[[ta_col]]
    pr_full <- df[[pr_col]]
    ok <- is.finite(ta_full) & is.finite(pr_full)
    df_ok <- df[ok, , drop = FALSE]

    ta <- df_ok[[ta_col]]
    pr <- df_ok[[pr_col]]

    if (length(ta) < min_samples) {
      warning("Skipping group '", grp, "' (need at least ", min_samples, " complete observations).", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    if (stats::sd(ta) < .Machine$double.eps || stats::sd(pr) < .Machine$double.eps) {
      warning("Skipping group '", grp, "' (near-zero variance in temperature or precipitation).", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    # Observation weights (external or uniform), renormalized after filtering
    w_pack <- .scenwgt_get_obs_weights(df_ok, weights_col = weights_col, grp = grp, tag = "COPULA")
    w_raw <- w_pack$w_raw
    w_obs <- w_pack$w_obs

    # Group scaling context
    sc_ctx <- .scenwgt_group_scale_context(
      ta = ta, pr = pr,
      grid_ctx = grid_ctx,
      grp = grp,
      verbose = verbose,
      tag = "COPULA",
      area_weight = area_weight
    )

    ta_use <- .scenwgt_apply_scale(ta, sc_ctx$scale_ta)
    pr_use <- .scenwgt_apply_scale(pr, sc_ctx$scale_pr)

    bw_min_use <- .scenwgt_scale_bw(bw_min, sc_ctx$scale_ta, sc_ctx$scale_pr)
    bw_max_use <- .scenwgt_scale_bw(bw_max, sc_ctx$scale_ta, sc_ctx$scale_pr)

    # Bandwidth selection (marginal KDEs)
    bw_use <- bw
    bw_method_used <- "user"

    if (is.null(bw_use)) {

      bw_data_ta <- .scenwgt_pick_bw_1d(ta_use, bw_method = bw_method, bw_min = 0, bw_max = Inf)
      bw_data_pr <- .scenwgt_pick_bw_1d(pr_use, bw_method = bw_method, bw_min = 0, bw_max = Inf)

      if (!is.finite(bw_data_ta) || !is.finite(bw_data_pr)) {
        warning("Skipping group '", grp, "' (marginal bandwidth selection failed).", call. = FALSE)
        skipped <- c(skipped, grp)
        next
      }

      # Grid-step floor (same intent as KDE)
      bw_grid <- c(k[1] * sc_ctx$dT_use, k[2] * sc_ctx$dP_use)
      if (!is.finite(bw_grid[1])) bw_grid[1] <- 0
      if (!is.finite(bw_grid[2])) bw_grid[2] <- 0

      bw_use <- pmax(bw_grid, alpha * c(bw_data_ta, bw_data_pr), na.rm = TRUE)

      bw_use <- pmax(bw_use, bw_min_use)
      bw_use <- pmin(bw_use, bw_max_use)

      bw_method_used <- bw_method

      if (!.scenwgt_validate_bw_vec_1d(bw_use)) {
        warning("Skipping group '", grp, "' (invalid marginal bandwidths after clamping).", call. = FALSE)
        skipped <- c(skipped, grp)
        next
      }

    } else {
      bw_use <- .scenwgt_scale_bw(bw_use, sc_ctx$scale_ta, sc_ctx$scale_pr)
      if (!.scenwgt_validate_bw_vec_1d(bw_use)) {
        warning("Skipping group '", grp, "' (invalid user bandwidths after scaling).", call. = FALSE)
        skipped <- c(skipped, grp)
        next
      }
    }

    bw_ta <- bw_use[1]
    bw_pr <- bw_use[2]

    # Diagnostics: effective sample size
    eff_n <- .scenwgt_effective_n(w_raw)

    if (diagnostics) {
      diag_bw[[grp]] <- c(ta = bw_ta, pr = bw_pr)
      diag_eff_n[[grp]] <- eff_n
      diag_scaling[[grp]] <- list(ta = sc_ctx$scale_ta, pr = sc_ctx$scale_pr)
      diag_bw_method[[grp]] <- bw_method_used
    }

    if (isTRUE(verbose)) {
      bw_label <- if (scale == "none") "bw" else "bw_scaled"
      message(sprintf("[COPULA] Group=%s | %s_%s = (ta=%.3f, pr=%.3f)",
                      grp, bw_label, bw_method_used, bw_ta, bw_pr))
    }

    if (is.finite(eff_n) && eff_n < min_samples && isTRUE(verbose)) {
      warning("Group '", grp, "': effective sample size (", round(eff_n, 1),
              ") is below min_samples (", min_samples, ").", call. = FALSE)
    }

    # -------------------------------------------------------------------------
    # Fit Gaussian copula on pseudo-observations derived from KDE CDF marginals
    # -------------------------------------------------------------------------

    u_ta_obs <- .scenwgt_kde1_cdf(ta_use, ta_use, w_obs, bw = bw_ta, chunk_size = chunk_size)
    u_pr_obs <- .scenwgt_kde1_cdf(pr_use, pr_use, w_obs, bw = bw_pr, chunk_size = chunk_size)

    u_ta_obs <- .scenwgt_clamp_unit(u_ta_obs)
    u_pr_obs <- .scenwgt_clamp_unit(u_pr_obs)

    z_ta <- stats::qnorm(u_ta_obs)
    z_pr <- stats::qnorm(u_pr_obs)

    rho <- .scenwgt_wcor(z_ta, z_pr, w_obs)

    if (!is.finite(rho)) {
      warning("Skipping group '", grp, "' (copula dependence fit failed: non-finite rho).", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    if (diagnostics) {
      diag_rho[[grp]] <- rho
    }

    if (isTRUE(verbose)) {
      message(sprintf("[COPULA] Group=%s | rho=%.3f", grp, rho))
    }

    # -------------------------------------------------------------------------
    # Evaluate joint density on the evaluation surface
    # -------------------------------------------------------------------------

    weights <- rep(NA_real_, n_grid)

    for (s in seq.int(1L, n_grid, by = chunk_size)) {
      e <- as.integer(min(s + chunk_size - 1L, n_grid))

      g_ta <- sc_ctx$grid_ta_use[s:e]
      g_pr <- sc_ctx$grid_pr_use[s:e]

      # Marginal pdfs and cdfs on the grid chunk
      f_ta <- .scenwgt_kde1_pdf(g_ta, ta_use, w_obs, bw = bw_ta, chunk_size = as.integer(e - s + 1L))
      f_pr <- .scenwgt_kde1_pdf(g_pr, pr_use, w_obs, bw = bw_pr, chunk_size = as.integer(e - s + 1L))

      u_ta <- .scenwgt_kde1_cdf(g_ta, ta_use, w_obs, bw = bw_ta, chunk_size = as.integer(e - s + 1L))
      u_pr <- .scenwgt_kde1_cdf(g_pr, pr_use, w_obs, bw = bw_pr, chunk_size = as.integer(e - s + 1L))

      u_ta <- .scenwgt_clamp_unit(u_ta)
      u_pr <- .scenwgt_clamp_unit(u_pr)

      z1 <- stats::qnorm(u_ta)
      z2 <- stats::qnorm(u_pr)

      c_dens <- .scenwgt_gaussian_copula_density(z1, z2, rho)

      weights[s:e] <- c_dens * f_ta * f_pr
    }

    if (normalize && sc_ctx$area_weight_use == "regular") {
      weights <- weights * (sc_ctx$dT_use * sc_ctx$dP_use)
    }

    # Support masking
    if (!is.null(support_mask)) {
      weights[!support_mask] <- 0
      if (all(weights == 0 | !is.finite(weights))) {
        warning("Group '", grp, "': all surface points outside support; weights are all zero.", call. = FALSE)
      }
    }

    if (normalize) weights <- .scenwgt_normalize_vec(weights)

    # Output col
    col_name <- .scenwgt_unique_colname(grp, existing_cols)
    existing_cols <- c(existing_cols, col_name)
    out[[col_name]] <- weights
  }

  out <- .scenwgt_order_weight_cols(out, scenario_grid)

  attr(out, "skipped_groups") <- if (length(skipped) > 0L) skipped else character(0)

  if (diagnostics) {
    attr(out, "bandwidth_used") <- diag_bw
    attr(out, "bandwidth_method") <- diag_bw_method
    attr(out, "copula_rho") <- diag_rho
    attr(out, "effective_sample_size") <- diag_eff_n
    attr(out, "scaling_params") <- list(
      mode = scale,
      global_ta = grid_ctx$scale_global_ta,
      global_pr = grid_ctx$scale_global_pr,
      by_group = diag_scaling
    )
  }

  out
}
