# ------------------------------------------------------------------------------
# KDE-specific helpers
# ------------------------------------------------------------------------------

#' Validate a 2D bandwidth vector
#'
#' @description
#' Internal predicate to validate a 2D KDE bandwidth vector `c(bw_ta, bw_pr)`.
#'
#' @details
#' A valid bandwidth vector is numeric, length 2, finite, and strictly positive.
#' Bandwidth units must match the units of the (possibly scaled) temperature and
#' precipitation inputs used by the KDE.
#'
#' @param bw_vec Numeric vector of length 2. Units: `[ta units, pr units]`.
#'
#' @return Logical scalar; `TRUE` if `bw_vec` is valid, otherwise `FALSE`.
#'
#' @keywords internal
.scenwgt_validate_bw_vec <- function(bw_vec) {
  is.numeric(bw_vec) && length(bw_vec) == 2L &&
    all(is.finite(bw_vec)) && all(bw_vec > 0)
}

#' Try 2D plug-in bandwidth via ks::Hpi
#'
#' @description
#' Attempts to estimate a 2D bandwidth matrix using a plug-in rule (via `ks::Hpi`)
#' and converts it to a length-2 bandwidth vector using `sqrt(diag(H))`.
#'
#' @details
#' - Input is treated as a bivariate sample `(ta, pr)` with observations in rows.
#' - Returns `NULL` if the plug-in estimator fails (e.g., due to numerical issues,
#'   degenerate covariance, or package-level errors).
#' - The returned bandwidths correspond to marginal smoothing scales implied by
#'   the diagonal of the 2D bandwidth matrix.
#'
#' @param ta Numeric vector of temperature values (finite). Units: `ta units`.
#' @param pr Numeric vector of precipitation values (finite). Units: `pr units`.
#'
#' @return Numeric vector `c(bw_ta, bw_pr)` or `NULL` if estimation fails.
#'
#' @keywords internal
.scenwgt_try_bw_plugin <- function(ta, pr) {
  data_mat <- cbind(ta, pr)
  H <- .kde_try_Hpi_2d(data_mat)
  if (is.null(H)) return(NULL)
  sqrt(diag(H))
}

#' Try 1D NRD bandwidth (diagonal 2D) for temperature and precipitation
#'
#' @description
#' Estimates separate 1D NRD bandwidths for `ta` and `pr` and returns them as a
#' diagonal 2D bandwidth vector `c(bw_ta, bw_pr)`.
#'
#' @details
#' This is a pragmatic fallback when plug-in multivariate bandwidth estimation is
#' unavailable or unstable. The underlying helper `.kde_try_nrd_2d()` is expected
#' to return a 2x2 diagonal bandwidth matrix or `NULL`.
#'
#' @param ta Numeric vector of temperature values (finite). Units: `ta units`.
#' @param pr Numeric vector of precipitation values (finite). Units: `pr units`.
#'
#' @return Numeric vector `c(bw_ta, bw_pr)` or `NULL` if estimation fails.
#'
#' @keywords internal
.scenwgt_try_bw_nrd <- function(ta, pr) {
  data_mat <- cbind(ta, pr)
  H <- .kde_try_nrd_2d(data_mat[, 1], data_mat[, 2])
  if (is.null(H)) return(NULL)
  sqrt(diag(H))
}

#' Select an automatic bandwidth with optional plug-in method and safety bounds
#'
#' @description
#' Chooses a bandwidth vector for 2D Gaussian product-kernel KDE by combining
#' (i) a grid-step-based floor and (ii) a data-driven estimate, then clamping to
#' user-specified min/max bounds.
#'
#' @details
#' The bandwidth is computed in three stages:
#' 1. **Grid floor:** `bw_grid = c(k[1] * dT, k[2] * dP)` where `dT`/`dP` are the
#'    grid steps (in the **current scale** used for evaluation).
#' 2. **Data-driven:** either plug-in (`ks::Hpi`) or NRD (per-dimension), selected
#'    by `bw_method`, with controlled fallback behavior when plug-in fails.
#' 3. **Combine + clamp:** `bw_use = pmax(bw_grid, alpha * bw_data)` then clamped
#'    elementwise to `[bw_min, bw_max]`.
#'
#' Returns `NULL` if a valid bandwidth cannot be produced after clamping/validation.
#'
#' Assumptions:
#' - `ta` and `pr` are finite and have non-negligible variance.
#' - `dT` and `dP` are non-negative and represent grid resolution in the same scale
#'   as `ta`/`pr` passed to this function.
#'
#' @param ta Numeric vector (finite) of temperature values in the evaluation scale.
#' @param pr Numeric vector (finite) of precipitation values in the evaluation scale.
#' @param dT Numeric scalar; grid step in the temperature dimension (>= 0). Units: `ta units`.
#' @param dP Numeric scalar; grid step in the precipitation dimension (>= 0). Units: `pr units`.
#' @param k Numeric length-2; multipliers for grid-step floor. Dimensionless.
#' @param alpha Numeric scalar > 0; multiplier applied to the data-driven bandwidth.
#' @param bw_min Numeric length-2; elementwise lower bounds (>= 0).
#' @param bw_max Numeric length-2; elementwise upper bounds (> 0 or `Inf`).
#' @param bw_method Character; one of `"auto"`, `"plugin"`, `"nrd"`.
#'
#' @return Numeric length-2 bandwidth vector `c(bw_ta, bw_pr)` or `NULL`.
#'
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
#'
#' @description
#' Computes weighted KDE values on a chunk of grid points using a separable
#' Gaussian kernel: `K(u, v) = phi(u) * phi(v)` with `u = (g_ta - ta)/bw_ta`,
#' `v = (g_pr - pr)/bw_pr`. This helper returns the unnormalized density sum
#' (i.e., it does **not** apply the `1/(bw_ta*bw_pr)` normalization constant).
#'
#' @details
#' Inputs must already be filtered to finite observations and consistent lengths.
#' Weights `w_obs` are expected to be normalized to sum to 1 (upstream logic does
#' this), but the computation itself works for any non-negative weights.
#'
#' Complexity:
#' - Time: `O(n_grid_chunk * n_obs)`
#' - Space: `O(n_grid_chunk * n_obs)` due to `outer()` allocations.
#'
#' @param g_ta Numeric vector of grid temperature coordinates for this chunk.
#'   Units: `ta units` (evaluation scale).
#' @param g_pr Numeric vector of grid precipitation coordinates for this chunk.
#'   Units: `pr units` (evaluation scale).
#' @param ta_use Numeric vector of observed temperatures (evaluation scale).
#' @param pr_use Numeric vector of observed precipitation values (evaluation scale).
#' @param w_obs Numeric vector of observation weights (same length as `ta_use`).
#' @param inv_bw_ta Numeric scalar; `1 / bw_ta`. Units: `1 / ta units`.
#' @param inv_bw_pr Numeric scalar; `1 / bw_pr`. Units: `1 / pr units`.
#'
#' @return Numeric vector of length `length(g_ta)` containing weighted kernel sums
#'   (without the `1/(bw_ta*bw_pr)` factor).
#'
#' @keywords internal
.scenwgt_kde_chunk_vectorized <- function(g_ta, g_pr, ta_use, pr_use,
                                          w_obs, inv_bw_ta, inv_bw_pr) {
  u_mat <- outer(g_ta, ta_use, `-`) * inv_bw_ta
  v_mat <- outer(g_pr, pr_use, `-`) * inv_bw_pr
  kernel_vals <- stats::dnorm(u_mat) * stats::dnorm(v_mat)
  as.numeric(kernel_vals %*% w_obs)
}

# ==============================================================================
# Main KDE function
# ==============================================================================

#' Calculate climate-informed scenario weights using 2D KDE
#'
#' @description
#' For each group (e.g., SSP), fits a 2D Gaussian product-kernel KDE to ensemble
#' points `(ta, pr)` and evaluates it on `scenario_grid`. If `normalize = TRUE`,
#' returns a discrete probability mass function (PMF) over grid points (optionally
#' scaled by grid-cell area for regular grids).
#'
#' @details
#' **What is computed**
#' - Let observations be `(ta_i, pr_i)` with weights `w_i` (either uniform or taken
#'   from `weights_col`, then renormalized to sum to 1 after filtering).
#' - For each grid point `(T, P)`, compute
#'   \deqn{ \hat{f}(T, P) = \frac{1}{bw_{ta}\,bw_{pr}} \sum_i w_i \, \phi\!\left(\frac{T-ta_i}{bw_{ta}}\right)\,
#'   \phi\!\left(\frac{P-pr_i}{bw_{pr}}\right) }
#'   where `phi()` is the standard normal density.
#'
#' **Scaling**
#' - `scale = "none"` evaluates in the original units.
#' - `scale = "global"` applies one global scaling (derived in grid context).
#' - `scale = "by_group"` applies scaling per group.
#' Scaling affects bandwidth selection and evaluation because it changes units.
#'
#' **Area weighting (only when `normalize = TRUE`)**
#' - If `area_weight = "regular"`, weights are multiplied by `(dT * dP)` prior to
#'   normalization, where `dT` and `dP` are grid steps in the evaluation scale.
#' - If `area_weight = "none"`, normalization is over raw evaluated densities.
#'
#' **Support masking**
#' - If `support` is provided, grid points outside support are set to zero after
#'   evaluation. If all grid points are outside support, a warning is emitted and
#'   the returned weight column for that group will be all zeros (after masking).
#'
#' **Diagnostics**
#' When `diagnostics = TRUE`, attributes are attached:
#' - `"bandwidth_used"`: per-group bandwidths used.
#' - `"bandwidth_method"`: `"user"` or the selected automatic method.
#' - `"effective_sample_size"`: Kish effective sample size computed from raw weights.
#' - `"scaling_params"`: scaling mode and per-group scaling factors.
#'
#' **Failure / skip conditions**
#' Groups are skipped (with a warning) when:
#' - fewer than `min_samples` complete `(ta, pr)` pairs remain,
#' - near-zero variance in `ta` or `pr`,
#' - automatic bandwidth selection fails.
#' Skipped groups are recorded in `attr(out, "skipped_groups")`.
#'
#' @param ensemble_data Data frame with ensemble projections. Must contain columns
#'   `ta_col`, `pr_col`, and `group_col`. Rows with non-finite `ta` or `pr` are ignored.
#' @param scenario_grid Data frame defining the evaluation grid. Must contain columns
#'   `ta_col` and `pr_col`. All rows are evaluated.
#' @param pr_col Character scalar; name of precipitation column in both inputs. Units: `pr units`.
#' @param ta_col Character scalar; name of temperature column in both inputs. Units: `ta units`.
#' @param group_col Character scalar; name of grouping column in `ensemble_data`
#'   (e.g., scenario/SSP). One output weight column is produced per group.
#' @param bw Optional numeric length-2 bandwidth vector `c(bw_ta, bw_pr)` in the
#'   **original units**. If `scale != "none"`, the bandwidth is internally scaled
#'   to the evaluation scale. Must be strictly positive.
#' @param bw_method Bandwidth selection method used when `bw` is `NULL`.
#'   One of `"auto"`, `"plugin"`, `"nrd"`. See Details.
#' @param k Numeric length-2; grid-step multipliers used for bandwidth floor.
#'   Dimensionless. Larger values enforce smoother KDE relative to grid resolution.
#' @param alpha Numeric scalar > 0; multiplier applied to the data-driven bandwidth
#'   estimate before combining with the grid floor.
#' @param bw_min Numeric length-2; elementwise lower bounds for bandwidth (>= 0) in
#'   the **original units**, later scaled if applicable.
#' @param bw_max Numeric length-2; elementwise upper bounds for bandwidth (> 0 or `Inf`)
#'   in the **original units**, later scaled if applicable.
#' @param min_samples Integer; minimum number of complete observations required per
#'   group after filtering. Must be >= 3.
#' @param normalize Logical; if `TRUE`, normalize each group’s weights to sum to 1
#'   (after optional area-weighting and support masking). If `FALSE`, returns raw
#'   evaluated densities (scaled by `1/(bw_ta*bw_pr)`).
#' @param area_weight Character; `"regular"` or `"none"`. Only applied when
#'   `normalize = TRUE`. See Details.
#' @param scale Character; `"none"`, `"global"`, or `"by_group"`. Controls scaling
#'   applied to the `(ta, pr)` space and grid before KDE evaluation.
#' @param weights_col Optional character scalar naming a column in `ensemble_data`
#'   containing per-row weights. Weights are validated and renormalized within each
#'   group after filtering to complete cases.
#' @param support Optional list specifying support bounds. Expected keys are `ta`
#'   and/or `pr`, each a numeric length-2 range `c(min, max)` in the **original units**.
#'   Grid points outside support are set to zero.
#' @param chunk_size Integer; number of grid rows evaluated per chunk. Larger values
#'   are faster but use more memory due to `outer()` allocations.
#' @param diagnostics Logical; if `TRUE`, attaches diagnostic attributes (see Details).
#' @param verbose Logical; if `TRUE`, prints per-group bandwidth info and emits
#'   additional warnings (e.g., low effective sample size).
#'
#' @return A data frame with the original `scenario_grid` columns plus one numeric
#'   weight column per group. Weight columns are named using the group labels,
#'   with collision avoidance if necessary.
#'
#' @examples
#' set.seed(1)
#' ensemble_data <- data.frame(
#'   scenario = rep(c("SSP1", "SSP2"), each = 30),
#'   tavg = c(rnorm(30, 0, 1), rnorm(30, 1, 1)),
#'   prcp = c(rnorm(30, 0, 1), rnorm(30, 0.5, 1))
#' )
#'
#' scenario_grid <- expand.grid(
#'   tavg = seq(-3, 3, length.out = 50),
#'   prcp = seq(-3, 3, length.out = 50)
#' )
#'
#' w <- compute_scenario_surface_weights_kde(
#'   ensemble_data = ensemble_data,
#'   scenario_grid = scenario_grid,
#'   ta_col = "tavg",
#'   pr_col = "prcp",
#'   group_col = "scenario",
#'   normalize = TRUE,
#'   area_weight = "none",
#'   verbose = FALSE
#' )
#' head(w)
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

  # Collect per-scenario long outputs and bind at the end.
  out_list <- list()
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

    # Long-format output: one row per grid point per scenario
    out_list[[length(out_list) + 1L]] <- .scenwgt_make_long_surface(
      scenario_grid = scenario_grid,
      scenario = as.character(grp),
      weight = weights
    )
  }

  if (length(out_list) == 0L) {
    out <- scenario_grid[0, , drop = FALSE]
    out$scenario <- character(0)
    out$weight <- numeric(0)
  } else {
    out <- do.call(rbind, out_list)
    rownames(out) <- NULL
  }

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

