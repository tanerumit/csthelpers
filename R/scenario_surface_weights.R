# ==============================================================================
# Climate Scenario Probability Estimation (REVISED)
# ==============================================================================

.scenwgt_grid_step <- function(x) {
  x <- sort(unique(x[is.finite(x)]))
  if (length(x) < 2L) return(NA_real_)
  stats::median(diff(x))
}

#' Normalize vector to sum to 1 (guarded)
#' @keywords internal
.scenwgt_normalize_vec <- function(x) {
  s <- sum(x)
  if (is.finite(s) && s > .Machine$double.eps) x / s else x
}

#' Compute scaling parameters (mean/sd) with near-zero sd guard
#' @keywords internal
.scenwgt_scale_params <- function(x, label, context, verbose, tag) {
  x <- x[is.finite(x)]
  mu <- mean(x)
  sd <- stats::sd(x)

  if (!is.finite(sd) || sd <= .Machine$double.eps) {
    if (isTRUE(verbose)) {
      message("[", tag, "] scaling disabled for ", label, context,
              " (non-finite or near-zero sd).")
    }
    return(list(mean = 0, sd = 1, enabled = FALSE))
  }

  list(mean = mu, sd = sd, enabled = TRUE)
}

#' Apply scaling if enabled
#' @keywords internal
.scenwgt_apply_scale <- function(x, params) {
  if (isTRUE(params$enabled)) (x - params$mean) / params$sd else x
}

#' Transform a step size into scaled space
#' @keywords internal
.scenwgt_scale_step <- function(step, params) {
  if (isTRUE(params$enabled)) step / params$sd else step
}

#' Convert bandwidths in original units to scaled units (per dimension)
#' @keywords internal
.scenwgt_scale_bw <- function(bw_vec, params_ta, params_pr) {
  bw_scaled <- bw_vec
  if (isTRUE(params_ta$enabled)) bw_scaled[1] <- bw_scaled[1] / params_ta$sd
  if (isTRUE(params_pr$enabled)) bw_scaled[2] <- bw_scaled[2] / params_pr$sd
  bw_scaled
}

#' Validate `support` structure (optional hard bounds on grid)
#' @keywords internal
.scenwgt_validate_support <- function(support) {
  if (is.null(support)) return(NULL)
  if (!is.list(support)) {
    stop("support must be a named list with elements 'ta' and/or 'pr'.", call. = FALSE)
  }
  if (length(support) == 0L) return(NULL)
  if (is.null(names(support)) || any(!nzchar(names(support)))) {
    stop("support must be a named list with elements 'ta' and/or 'pr'.", call. = FALSE)
  }

  bad_names <- setdiff(names(support), c("ta", "pr"))
  if (length(bad_names) > 0L) {
    stop("support can only contain 'ta' and/or 'pr'.", call. = FALSE)
  }

  .check_range <- function(x, label) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
      stop("support$", label, " must be numeric length-2 with finite values.", call. = FALSE)
    }
    if (x[1] > x[2]) stop("support$", label, " must have min <= max.", call. = FALSE)
    x
  }

  if (!is.null(support$ta)) support$ta <- .check_range(support$ta, "ta")
  if (!is.null(support$pr)) support$pr <- .check_range(support$pr, "pr")
  support
}

#' Boolean mask of grid points inside `support`
#' @keywords internal
.scenwgt_support_mask <- function(grid_ta, grid_pr, support) {
  if (is.null(support)) return(NULL)
  mask <- rep(TRUE, length(grid_ta))

  if (!is.null(support$ta)) {
    mask <- mask & is.finite(grid_ta) & grid_ta >= support$ta[1] & grid_ta <= support$ta[2]
  }
  if (!is.null(support$pr)) {
    mask <- mask & is.finite(grid_pr) & grid_pr >= support$pr[1] & grid_pr <= support$pr[2]
  }

  mask
}

#' Disable regular-grid area weighting if grid steps are invalid
#' @keywords internal
.scenwgt_area_weight_guard <- function(area_weight, dT, dP, verbose, tag) {
  if (area_weight != "regular") return(area_weight)

  if (!is.finite(dT) || !is.finite(dP) ||
      dT <= .Machine$double.eps || dP <= .Machine$double.eps) {
    if (isTRUE(verbose)) message("[", tag, "] area_weight='regular' disabled (invalid grid step).")
    return("none")
  }

  "regular"
}

#' Prepare global grid context: steps, support mask, and optional global scaling
#' @keywords internal
.scenwgt_prepare_grid_context <- function(
    ensemble_data, scenario_grid,
    ta_col, pr_col,
    scale, area_weight, support,
    verbose, tag) {

  scale <- match.arg(scale, c("none", "global", "by_group"))
  area_weight <- match.arg(area_weight, c("regular", "none"))
  support <- .scenwgt_validate_support(support)

  grid_ta <- scenario_grid[[ta_col]]
  grid_pr <- scenario_grid[[pr_col]]
  n_grid <- nrow(scenario_grid)

  dT <- .scenwgt_grid_step(grid_ta)
  dP <- .scenwgt_grid_step(grid_pr)

  support_mask <- .scenwgt_support_mask(grid_ta, grid_pr, support)

  scale_none <- list(mean = 0, sd = 1, enabled = FALSE)

  scale_global_ta <- scale_none
  scale_global_pr <- scale_none

  grid_ta_global <- grid_ta
  grid_pr_global <- grid_pr
  dT_global <- dT
  dP_global <- dP

  if (scale == "global") {
    scale_global_ta <- .scenwgt_scale_params(ensemble_data[[ta_col]], "ta", " (global)", verbose, tag)
    scale_global_pr <- .scenwgt_scale_params(ensemble_data[[pr_col]], "pr", " (global)", verbose, tag)

    grid_ta_global <- .scenwgt_apply_scale(grid_ta, scale_global_ta)
    grid_pr_global <- .scenwgt_apply_scale(grid_pr, scale_global_pr)

    dT_global <- .scenwgt_scale_step(dT, scale_global_ta)
    dP_global <- .scenwgt_scale_step(dP, scale_global_pr)
  }

  area_weight_global <- .scenwgt_area_weight_guard(area_weight, dT_global, dP_global, verbose, tag)

  list(
    grid_ta = grid_ta,
    grid_pr = grid_pr,
    n_grid = n_grid,
    dT = dT,
    dP = dP,
    support = support,
    support_mask = support_mask,
    scale = scale,
    scale_none = scale_none,
    scale_global_ta = scale_global_ta,
    scale_global_pr = scale_global_pr,
    grid_ta_global = grid_ta_global,
    grid_pr_global = grid_pr_global,
    dT_global = dT_global,
    dP_global = dP_global,
    area_weight_global = area_weight_global
  )
}

#' Prepare per-group scaling view of the grid and steps
#' @keywords internal
.scenwgt_group_scale_context <- function(
    ta, pr, grid_ctx, grp, verbose, tag, area_weight) {

  scale <- grid_ctx$scale

  if (scale == "none") {
    scale_ta <- grid_ctx$scale_none
    scale_pr <- grid_ctx$scale_none
    grid_ta_use <- grid_ctx$grid_ta
    grid_pr_use <- grid_ctx$grid_pr
    dT_use <- grid_ctx$dT
    dP_use <- grid_ctx$dP
    area_weight_use <- grid_ctx$area_weight_global

  } else if (scale == "global") {
    scale_ta <- grid_ctx$scale_global_ta
    scale_pr <- grid_ctx$scale_global_pr
    grid_ta_use <- grid_ctx$grid_ta_global
    grid_pr_use <- grid_ctx$grid_pr_global
    dT_use <- grid_ctx$dT_global
    dP_use <- grid_ctx$dP_global
    area_weight_use <- grid_ctx$area_weight_global

  } else {
    scale_ta <- .scenwgt_scale_params(ta, "ta", paste0(" (group=", grp, ")"), verbose, tag)
    scale_pr <- .scenwgt_scale_params(pr, "pr", paste0(" (group=", grp, ")"), verbose, tag)

    grid_ta_use <- .scenwgt_apply_scale(grid_ctx$grid_ta, scale_ta)
    grid_pr_use <- .scenwgt_apply_scale(grid_ctx$grid_pr, scale_pr)

    dT_use <- .scenwgt_scale_step(grid_ctx$dT, scale_ta)
    dP_use <- .scenwgt_scale_step(grid_ctx$dP, scale_pr)

    area_weight_use <- .scenwgt_area_weight_guard(area_weight, dT_use, dP_use, verbose, tag)
  }

  list(
    scale_ta = scale_ta,
    scale_pr = scale_pr,
    grid_ta_use = grid_ta_use,
    grid_pr_use = grid_pr_use,
    dT_use = dT_use,
    dP_use = dP_use,
    area_weight_use = area_weight_use
  )
}

#' Reorder output so base grid cols first, weight cols sorted
#' @keywords internal
.scenwgt_order_weight_cols <- function(out, scenario_grid) {
  base_cols <- names(scenario_grid)
  weight_cols <- setdiff(names(out), base_cols)
  out[, c(base_cols, sort(weight_cols)), drop = FALSE]
}

#' Extract and validate per-observation weights for a group (or uniform)
#' @keywords internal
.scenwgt_get_obs_weights <- function(df_ok, weights_col, grp, tag) {

  if (is.null(weights_col)) {
    n <- nrow(df_ok)
    w_raw <- rep(1, n)
    w_obs <- rep(1 / n, n)
    return(list(w_raw = w_raw, w_obs = w_obs))
  }

  w_raw <- df_ok[[weights_col]]

  if (!is.numeric(w_raw)) {
    stop("weights_col '", weights_col, "' must be numeric.", call. = FALSE)
  }
  if (length(w_raw) != nrow(df_ok)) {
    stop("Internal error: weights length mismatch after filtering in group '", grp, "'.", call. = FALSE)
  }
  if (any(!is.finite(w_raw)) || any(w_raw < 0)) {
    stop("weights_col '", weights_col, "' contains non-finite or negative values in group '", grp, "'.",
         call. = FALSE)
  }
  if (all(w_raw == 0)) {
    stop("weights_col '", weights_col, "' is all zero in group '", grp, "'.", call. = FALSE)
  }

  w_obs <- w_raw / sum(w_raw)

  list(w_raw = w_raw, w_obs = w_obs)
}

#' Effective sample size (Kish) for diagnostics
#' @keywords internal
.scenwgt_effective_n <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

#' Generate unique column name, checking for duplicates
#' @keywords internal
.scenwgt_unique_colname <- function(grp, existing_names) {
  col_name <- make.names(grp)
  if (col_name %in% existing_names) {
    stop("Duplicate column name '", col_name, "' after make.names() transformation ",
         "(from group '", grp, "'). Ensure group names produce unique column names.",
         call. = FALSE)
  }
  col_name
}


# ------------------------------------------------------------------------------
# KDE-specific helpers (REVISED)
# ------------------------------------------------------------------------------

#' Validate a 2D bandwidth vector
#' @keywords internal
.scenwgt_validate_bw_vec <- function(bw_vec) {
  is.numeric(bw_vec) && length(bw_vec) == 2L &&
    all(is.finite(bw_vec)) && all(bw_vec > 0)
}

#' Try 2D plug-in bandwidth via ks::Hpi
#' @keywords internal
.scenwgt_try_bw_plugin <- function(ta, pr) {
  data_mat <- cbind(ta, pr)
  H <- .kde_try_Hpi_2d(data_mat)
  if (is.null(H)) return(NULL)
  sqrt(diag(H))
}


#' Try 1D NRD bandwidth
#' @keywords internal
.scenwgt_try_bw_nrd <- function(ta, pr) {
  data_mat <- cbind(ta, pr)
  H <- .kde_try_nrd_2d(data_mat[, 1], data_mat[, 2])
  if (is.null(H)) return(NULL)
  sqrt(diag(H))
}


#' Auto bandwidth with optional plug-in method
#' @keywords internal
.scenwgt_pick_bw_auto <- function(ta, pr, dT, dP, k, alpha, bw_min, bw_max,
                                   bw_method = c("auto", "plugin", "nrd")) {
  bw_method <- match.arg(bw_method)

  # Grid-step floor
  bw_grid <- c(k[1] * dT, k[2] * dP)
  if (!is.finite(bw_grid[1])) bw_grid[1] <- 0
  if (!is.finite(bw_grid[2])) bw_grid[2] <- 0

  # Data-driven bandwidth
  bw_data <- NULL

  if (bw_method == "plugin") {
    bw_data <- .scenwgt_try_bw_plugin(ta, pr)
    if (is.null(bw_data)) {
      warning("ks::Hpi failed; falling back to NRD bandwidth.", call. = FALSE)
      bw_data <- .scenwgt_try_bw_nrd(ta, pr)
    }
  } else if (bw_method == "nrd") {
    bw_data <- .scenwgt_try_bw_nrd(ta, pr)
  } else {
    # auto: try plugin first, fall back to nrd
    bw_data <- .scenwgt_try_bw_plugin(ta, pr)
    if (is.null(bw_data)) {
      bw_data <- .scenwgt_try_bw_nrd(ta, pr)
    }
  }

  if (is.null(bw_data)) {
    bw_data <- c(NA_real_, NA_real_)
  }

  # Combine: max(grid floor, alpha * data-driven)
  bw_use <- pmax(bw_grid, alpha * bw_data, na.rm = TRUE)

  # Clamp to [bw_min, bw_max]
  bw_use <- pmax(bw_use, bw_min)
  bw_use <- pmin(bw_use, bw_max)

  if (!.scenwgt_validate_bw_vec(bw_use)) return(NULL)
  bw_use
}

#' Vectorized weighted KDE over a chunk of grid points (Gaussian product kernel)
#' @keywords internal
.scenwgt_kde_chunk_vectorized <- function(g_ta, g_pr, ta_use, pr_use,
                                           w_obs, inv_bw_ta, inv_bw_pr) {
  u_mat <- outer(g_ta, ta_use, `-`) * inv_bw_ta
  v_mat <- outer(g_pr, pr_use, `-`) * inv_bw_pr
  kernel_vals <- stats::dnorm(u_mat) * stats::dnorm(v_mat)
  as.numeric(kernel_vals %*% w_obs)
}


# ------------------------------------------------------------------------------
# MVN-specific helpers (REVISED)
# ------------------------------------------------------------------------------

#' Weighted mean (guarded, accepts pre-normalized or raw weights)
#' @keywords internal
.scenwgt_wmean <- function(x, w) {
  s <- sum(w)
  if (!is.finite(s) || s <= .Machine$double.eps) return(NA_real_)
  sum(w * x) / s
}

#' Weighted 2x2 covariance for (x,y) with weights w
#' @keywords internal
.scenwgt_wcov2 <- function(x, y, w) {
  s <- sum(w)
  if (!is.finite(s) || s <= .Machine$double.eps) return(matrix(NA_real_, 2, 2))
  w <- w / s
  mx <- sum(w * x)
  my <- sum(w * y)
  dx <- x - mx
  dy <- y - my
  cxx <- sum(w * dx * dx)
  cyy <- sum(w * dy * dy)
  cxy <- sum(w * dx * dy)
  matrix(c(cxx, cxy, cxy, cyy), nrow = 2, byrow = TRUE)
}

#' Validate Sigma as finite, positive-definite, and well-conditioned
#' @keywords internal
.scenwgt_validate_sigma <- function(Sigma, grp = "") {
  if (any(!is.finite(Sigma))) {
    return(list(valid = FALSE, reason = "covariance matrix contains non-finite values",
                chol = NULL))
  }
  if (any(diag(Sigma) <= .Machine$double.eps)) {
    return(list(valid = FALSE, reason = "near-zero variance in fitted covariance",
                chol = NULL))
  }

  # FIX 4.1: Check condition number for near-singular matrices
  cond_num <- tryCatch(kappa(Sigma), error = function(e) Inf)
  if (!is.finite(cond_num) || cond_num > 1e10) {
    return(list(valid = FALSE, reason = "near-singular covariance (high condition number)",
                chol = NULL))
  }

  # Compute Cholesky and return it for reuse
  chol_result <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_result)) {
    return(list(valid = FALSE, reason = "singular or non-positive-definite covariance",
                chol = NULL))
  }

  list(valid = TRUE, reason = NULL, chol = chol_result)
}

#' Bivariate normal density using pre-computed Cholesky factor
#' @keywords internal
.scenwgt_dmvnorm2_chol <- function(x, mu, chol_L) {
  n <- nrow(x)

  log_det <- 2 * sum(log(diag(chol_L)))
  if (!is.finite(log_det)) return(rep(NA_real_, n))

  xc <- cbind(x[, 1] - mu[1], x[, 2] - mu[2])
  z <- backsolve(chol_L, t(xc), transpose = TRUE)
  qf <- colSums(z^2)

  log_dens <- -0.5 * (2 * log(2 * pi) + log_det + qf)
  exp(log_dens)
}


# ==============================================================================
# Main KDE function (REVISED)
# ==============================================================================

#' Calculate Climate-Informed Scenario Weights (KDE; external per-row weights)
#'
#' @description
#' For each group (e.g., SSP), fits a 2D Gaussian product-kernel KDE to ensemble
#' (ta, pr) points and evaluates it on `scenario_grid`. If `normalize=TRUE`, returns
#' a discrete PMF over grid points (optionally area-weighted for regular grids).
#'
#' @param ensemble_data Data frame containing ensemble projections.
#' @param scenario_grid Data frame defining the evaluation grid.
#' @param pr_col Column name for precipitation variable.
#' @param ta_col Column name for temperature variable.
#' @param group_col Column name for grouping (e.g., scenario).
#' @param bw Optional user-specified bandwidth vector c(bw_ta, bw_pr).
#' @param bw_method Bandwidth selection method: "auto", "plugin", or "nrd".
#' @param k Grid-step multipliers for bandwidth floor.
#' @param alpha Multiplier for data-driven bandwidth.
#' @param bw_min Minimum bandwidth bounds.
#' @param bw_max Maximum bandwidth bounds.
#' @param min_samples Minimum observations required per group.
#' @param normalize Logical; normalize to PMF?
#' @param area_weight Area weighting method: "regular" or "none".
#' @param scale Scaling method: "none", "global", or "by_group".
#' @param weights_col Optional column with per-row weights.
#' @param support Optional list with ta and/or pr bounds.
#' @param chunk_size Grid evaluation chunk size.
#' @param diagnostics Logical; attach diagnostic attributes?
#' @param verbose Logical; print progress messages?
#'
#' @return Data frame with scenario_grid columns plus one weight column per group.
#'
#' @importFrom MASS bandwidth.nrd
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
compute_scenario_surface_weights_kde <- function(
    ensemble_data,
    scenario_grid,
    pr_col = "prcp",
    ta_col = "tavg",
    group_col = "scenario",
    bw = NULL,
    bw_method = c("auto", "plugin", "nrd"),
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

  if (!is.logical(normalize) || length(normalize) != 1L) stop("normalize must be a single logical.", call. = FALSE)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L) stop("diagnostics must be a single logical.", call. = FALSE)

  area_weight <- match.arg(area_weight)
  scale <- match.arg(scale)
  bw_method <- match.arg(bw_method)

  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 3L) {
    stop("min_samples must be a single integer >= 3.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  if (!is.null(bw) && !.scenwgt_validate_bw_vec(bw)) {
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
    tag = "KDE"
  )

  n_grid <- grid_ctx$n_grid
  support_mask <- grid_ctx$support_mask

  # ---------------------------------------------------------------------------
  # Loop over groups (FIX 3.2: use split())
  # ---------------------------------------------------------------------------

  out <- scenario_grid
  existing_cols <- names(out)
  skipped <- character(0)

  diag_bandwidth <- list()
  diag_eff_n <- list()
  diag_scaling <- list()
  diag_bw_method <- list()

  # Split data by group for efficiency
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
    w_pack <- .scenwgt_get_obs_weights(df_ok, weights_col = weights_col, grp = grp, tag = "KDE")
    w_raw <- w_pack$w_raw
    w_obs <- w_pack$w_obs

    # Group scaling context
    sc_ctx <- .scenwgt_group_scale_context(
      ta = ta, pr = pr,
      grid_ctx = grid_ctx,
      grp = grp,
      verbose = verbose,
      tag = "KDE",
      area_weight = area_weight
    )

    ta_use <- .scenwgt_apply_scale(ta, sc_ctx$scale_ta)
    pr_use <- .scenwgt_apply_scale(pr, sc_ctx$scale_pr)

    bw_min_use <- .scenwgt_scale_bw(bw_min, sc_ctx$scale_ta, sc_ctx$scale_pr)
    bw_max_use <- .scenwgt_scale_bw(bw_max, sc_ctx$scale_ta, sc_ctx$scale_pr)

    # Bandwidth selection
    bw_use <- bw
    bw_method_used <- "user"

    if (is.null(bw_use)) {
      bw_use <- .scenwgt_pick_bw_auto(
        ta = ta_use, pr = pr_use,
        dT = sc_ctx$dT_use, dP = sc_ctx$dP_use,
        k = k, alpha = alpha,
        bw_min = bw_min_use, bw_max = bw_max_use,
        bw_method = bw_method
      )
      bw_method_used <- bw_method

      if (is.null(bw_use)) {
        warning("Skipping group '", grp, "' (auto bandwidth selection failed).", call. = FALSE)
        skipped <- c(skipped, grp)
        next
      }
    } else {
      bw_use <- .scenwgt_scale_bw(bw_use, sc_ctx$scale_ta, sc_ctx$scale_pr)
    }

    if (isTRUE(verbose)) {
      bw_label <- if (scale == "none") "bw" else "bw_scaled"
      message(sprintf("[KDE] Group=%s | %s_%s = (ta=%.3f, pr=%.3f)",
                      grp, bw_label, bw_method_used, bw_use[1], bw_use[2]))
    }

    bw_ta <- bw_use[1]
    bw_pr <- bw_use[2]

    # Diagnostics: check effective sample size
    eff_n <- .scenwgt_effective_n(w_raw)
    if (diagnostics) {
      diag_bandwidth[[grp]] <- c(ta = bw_ta, pr = bw_pr)
      diag_eff_n[[grp]] <- eff_n
      diag_scaling[[grp]] <- list(ta = sc_ctx$scale_ta, pr = sc_ctx$scale_pr)
      diag_bw_method[[grp]] <- bw_method_used
    }

    # Warn if effective N is very low
    if (is.finite(eff_n) && eff_n < min_samples && isTRUE(verbose)) {
      warning("Group '", grp, "': effective sample size (", round(eff_n, 1),
              ") is below min_samples (", min_samples, ").", call. = FALSE)
    }

    # KDE evaluation (chunked)
    dens <- numeric(n_grid)
    inv_bw_ta <- 1 / bw_ta
    inv_bw_pr <- 1 / bw_pr
    norm_const <- inv_bw_ta * inv_bw_pr

    for (s in seq.int(1L, n_grid, by = chunk_size)) {
      e <- as.integer(min(s + chunk_size - 1L, n_grid))
      g_ta <- sc_ctx$grid_ta_use[s:e]
      g_pr <- sc_ctx$grid_pr_use[s:e]

      chunk_dens <- .scenwgt_kde_chunk_vectorized(
        g_ta = g_ta,
        g_pr = g_pr,
        ta_use = ta_use,
        pr_use = pr_use,
        w_obs = w_obs,
        inv_bw_ta = inv_bw_ta,
        inv_bw_pr = inv_bw_pr
      )

      dens[s:e] <- chunk_dens * norm_const
    }

    weights <- dens
    if (normalize && sc_ctx$area_weight_use == "regular") {
      weights <- weights * (sc_ctx$dT_use * sc_ctx$dP_use)
    }

    # FIX 4.3: Warn if all weights are zero after support masking
    if (!is.null(support_mask)) {
      weights[!support_mask] <- 0
      if (all(weights == 0 | !is.finite(weights))) {
        warning("Group '", grp, "': all grid points outside support; weights are all zero.", call. = FALSE)
      }
    }

    if (normalize) weights <- .scenwgt_normalize_vec(weights)

    # FIX 2.2: Check for duplicate column names
    col_name <- .scenwgt_unique_colname(grp, existing_cols)
    existing_cols <- c(existing_cols, col_name)
    out[[col_name]] <- weights
  }

  out <- .scenwgt_order_weight_cols(out, scenario_grid)

  # Attributes
  attr(out, "skipped_groups") <- if (length(skipped) > 0L) skipped else character(0)

  if (diagnostics) {
    attr(out, "bandwidth_used") <- diag_bandwidth
    attr(out, "bandwidth_method") <- diag_bw_method
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


# ==============================================================================
# Main MVN function (REVISED)
# ==============================================================================

#' Calculate Climate-Informed Scenario Weights (MVN; external per-row weights)
#'
#' @description
#' For each group, fits a bivariate normal distribution to (ta, pr) and evaluates
#' its density on `scenario_grid`. If `normalize=TRUE`, returns a discrete PMF over
#' grid points (optionally area-weighted for regular grids).
#'
#' @param ensemble_data Data frame containing ensemble projections.
#' @param scenario_grid Data frame defining the evaluation grid.
#' @param pr_col Column name for precipitation variable.
#' @param ta_col Column name for temperature variable.
#' @param group_col Column name for grouping (e.g., scenario).
#' @param robust Logical; use robust covariance estimation?
#' @param min_samples Minimum observations required per group.
#' @param normalize Logical; normalize to PMF?
#' @param area_weight Area weighting method: "regular" or "none".
#' @param scale Scaling method: "none", "global", or "by_group".
#' @param weights_col Optional column with per-row weights (incompatible with robust=TRUE).
#' @param support Optional list with ta and/or pr bounds.
#' @param diagnostics Logical; attach diagnostic attributes?
#' @param verbose Logical; print progress messages?
#'
#' @return Data frame with scenario_grid columns plus one weight column per group.
#'
#' @importFrom MASS cov.trob
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
compute_scenario_surface_weights_mvnorm <- function(
    ensemble_data,
    scenario_grid,
    pr_col = "prcp",
    ta_col = "tavg",
    group_col = "scenario",
    robust = TRUE,
    min_samples = 5L,
    normalize = TRUE,
    area_weight = c("regular", "none"),
    scale = c("none", "global", "by_group"),
    weights_col = NULL,
    support = NULL,
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

  if (!is.logical(robust) || length(robust) != 1L) stop("robust must be a single logical.", call. = FALSE)
  if (!is.logical(normalize) || length(normalize) != 1L) stop("normalize must be a single logical.", call. = FALSE)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L) stop("diagnostics must be a single logical.", call. = FALSE)

  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 3L) {
    stop("min_samples must be a single integer >= 3.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  area_weight <- match.arg(area_weight)
  scale <- match.arg(scale)

  if (!is.null(weights_col)) {
    if (!is.character(weights_col) || length(weights_col) != 1L || !nzchar(weights_col)) {
      stop("weights_col must be NULL or a single non-empty character column name.", call. = FALSE)
    }
    if (!(weights_col %in% names(ensemble_data))) {
      stop("weights_col not found in ensemble_data: ", weights_col, call. = FALSE)
    }
    if (isTRUE(robust)) {
      stop("robust=TRUE is not compatible with weights_col (cov.trob() does not support weights).",
           call. = FALSE)
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
    tag = "MVN"
  )

  n_grid <- grid_ctx$n_grid
  support_mask <- grid_ctx$support_mask

  # ---------------------------------------------------------------------------
  # Loop over groups (FIX 3.2: use split())
  # ---------------------------------------------------------------------------

  out <- scenario_grid
  existing_cols <- names(out)
  skipped <- character(0)

  diag_mvn_params <- list()
  diag_eff_n <- list()
  diag_scaling <- list()
  diag_fit_method <- list()

  # Split data by group for efficiency
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

    # Group scaling context
    sc_ctx <- .scenwgt_group_scale_context(
      ta = ta, pr = pr,
      grid_ctx = grid_ctx,
      grp = grp,
      verbose = verbose,
      tag = "MVN",
      area_weight = area_weight
    )

    ta_use <- .scenwgt_apply_scale(ta, sc_ctx$scale_ta)
    pr_use <- .scenwgt_apply_scale(pr, sc_ctx$scale_pr)

    # Observation weights (external or uniform), renormalized after filtering
    w_pack <- .scenwgt_get_obs_weights(df_ok, weights_col = weights_col, grp = grp, tag = "MVN")
    w_raw <- w_pack$w_raw
    w_obs <- w_pack$w_obs

    # Fit MVN
    fit <- tryCatch({
      if (isTRUE(robust)) {
        m <- MASS::cov.trob(cbind(ta_use, pr_use))
        list(mu = as.numeric(m$center), Sigma = m$cov, method = "robust")
      } else {
        mu <- c(.scenwgt_wmean(ta_use, w_obs), .scenwgt_wmean(pr_use, w_obs))
        Sigma <- .scenwgt_wcov2(ta_use, pr_use, w_obs)
        list(mu = mu, Sigma = Sigma, method = "classical_weighted")
      }
    }, error = function(e) {
      warning("Skipping group '", grp, "' (MVN fit failed: ", e$message, ").", call. = FALSE)
      NULL
    })

    if (is.null(fit)) {
      skipped <- c(skipped, grp)
      next
    }

    mu <- fit$mu
    Sigma <- fit$Sigma

    # FIX 2.1 & 3.4: Validate sigma with condition number check; reuse Cholesky
    sigma_check <- .scenwgt_validate_sigma(Sigma, grp)
    if (!sigma_check$valid) {
      warning("Skipping group '", grp, "' (", sigma_check$reason, ").", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    if (isTRUE(verbose)) {
      message(sprintf("[MVN] Group=%s | fit=%s | mu=(ta=%.3f, pr=%.3f)",
                      grp, fit$method, mu[1], mu[2]))
    }

    # Diagnostics: check effective sample size
    eff_n <- .scenwgt_effective_n(w_raw)
    if (diagnostics) {
      diag_mvn_params[[grp]] <- list(mu = mu, Sigma = Sigma)
      diag_fit_method[[grp]] <- fit$method
      diag_eff_n[[grp]] <- eff_n
      diag_scaling[[grp]] <- list(ta = sc_ctx$scale_ta, pr = sc_ctx$scale_pr)
    }

    # Warn if effective N is very low
    if (is.finite(eff_n) && eff_n < min_samples && isTRUE(verbose)) {
      warning("Group '", grp, "': effective sample size (", round(eff_n, 1),
              ") is below min_samples (", min_samples, ").", call. = FALSE)
    }

    # FIX 3.4: Reuse Cholesky factor from validation
    x_grid <- cbind(sc_ctx$grid_ta_use, sc_ctx$grid_pr_use)
    dens <- .scenwgt_dmvnorm2_chol(x_grid, mu = mu, chol_L = sigma_check$chol)

    if (all(!is.finite(dens))) {
      warning("Skipping group '", grp, "' (density evaluation failed).", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    weights <- dens
    if (normalize && sc_ctx$area_weight_use == "regular") {
      weights <- weights * (sc_ctx$dT_use * sc_ctx$dP_use)
    }

    # FIX 4.3: Warn if all weights are zero after support masking
    if (!is.null(support_mask)) {
      weights[!support_mask] <- 0
      if (all(weights == 0 | !is.finite(weights))) {
        warning("Group '", grp, "': all grid points outside support; weights are all zero.", call. = FALSE)
      }
    }

    if (normalize) weights <- .scenwgt_normalize_vec(weights)

    # FIX 2.2: Check for duplicate column names
    col_name <- .scenwgt_unique_colname(grp, existing_cols)
    existing_cols <- c(existing_cols, col_name)
    out[[col_name]] <- weights
  }

  out <- .scenwgt_order_weight_cols(out, scenario_grid)

  # Attributes
  attr(out, "skipped_groups") <- if (length(skipped) > 0L) skipped else character(0)

  if (diagnostics) {
    attr(out, "mvn_params") <- diag_mvn_params
    attr(out, "effective_sample_size") <- diag_eff_n
    attr(out, "scaling_params") <- list(
      mode = scale,
      global_ta = grid_ctx$scale_global_ta,
      global_pr = grid_ctx$scale_global_pr,
      by_group = diag_scaling
    )
    attr(out, "fit_method") <- diag_fit_method
  }

  out
}
