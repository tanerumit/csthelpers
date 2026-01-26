# ==============================================================================
# Climate Scenario Probability Estimation
# Revised with: Vectorized KDE, Cholesky-based MVN, Diagnostics Output
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared helper functions
# ------------------------------------------------------------------------------

.scenprob_grid_step <- function(x) {

  x <- sort(unique(x[is.finite(x)]))
  if (length(x) < 2) return(NA_real_)
  stats::median(diff(x))
}

.scenprob_normalize_vec <- function(x) {
  s <- sum(x)
  if (is.finite(s) && s > .Machine$double.eps) x / s else x
}

.scenprob_compute_family_weights <- function(family_id) {
  family_id <- as.character(family_id)
  family_id[is.na(family_id)] <- "NA_FAMILY"
  fam_tab <- table(family_id)
  n_fam <- length(fam_tab)
  w_fam <- 1 / n_fam
  w_i <- w_fam / as.numeric(fam_tab[family_id]) # split within family
  as.numeric(w_i)
}

.scenprob_scale_params <- function(x, label, context, verbose, tag) {
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

.scenprob_apply_scale <- function(x, params) {
  if (isTRUE(params$enabled)) (x - params$mean) / params$sd else x
}

.scenprob_scale_step <- function(step, params) {
  if (isTRUE(params$enabled)) step / params$sd else step
}

.scenprob_scale_bw <- function(bw_vec, params_ta, params_pr) {
  bw_scaled <- bw_vec
  if (isTRUE(params_ta$enabled)) bw_scaled[1] <- bw_scaled[1] / params_ta$sd
  if (isTRUE(params_pr$enabled)) bw_scaled[2] <- bw_scaled[2] / params_pr$sd
  bw_scaled
}

.scenprob_validate_support <- function(support) {
  if (is.null(support)) return(NULL)
  if (!is.list(support)) {
    stop("support must be a named list with elements 'ta' and/or 'pr'.", call. = FALSE)
  }
  if (length(support) == 0) return(NULL)
  if (is.null(names(support)) || any(!nzchar(names(support)))) {
    stop("support must be a named list with elements 'ta' and/or 'pr'.", call. = FALSE)
  }
  bad_names <- setdiff(names(support), c("ta", "pr"))
  if (length(bad_names) > 0) {
    stop("support can only contain 'ta' and/or 'pr'.", call. = FALSE)
  }

  .check_range <- function(x, label) {
    if (!is.numeric(x) || length(x) != 2 || any(!is.finite(x))) {
      stop("support$", label, " must be numeric length-2 with finite values.", call. = FALSE)
    }
    if (x[1] > x[2]) stop("support$", label, " must have min <= max.", call. = FALSE)
    x
  }

  if (!is.null(support$ta)) support$ta <- .check_range(support$ta, "ta")
  if (!is.null(support$pr)) support$pr <- .check_range(support$pr, "pr")
  support
}

.scenprob_support_mask <- function(grid_ta, grid_pr, support) {
  if (is.null(support)) return(NULL)
  n <- length(grid_ta)
  mask <- rep(TRUE, n)
  if (!is.null(support$ta)) {
    mask <- mask & is.finite(grid_ta) & grid_ta >= support$ta[1] & grid_ta <= support$ta[2]
  }
  if (!is.null(support$pr)) {
    mask <- mask & is.finite(grid_pr) & grid_pr >= support$pr[1] & grid_pr <= support$pr[2]
  }
  mask
}

.scenprob_area_weight_guard <- function(area_weight, dT, dP, verbose, tag) {
  if (area_weight != "regular") return(area_weight)
  if (!is.finite(dT) || !is.finite(dP) || dT <= .Machine$double.eps || dP <= .Machine$double.eps) {
    if (isTRUE(verbose)) message("[", tag, "] area_weight='regular' disabled (invalid grid step).")
    return("none")
  }
  "regular"
}

.scenprob_prepare_grid_context <- function(
    ensemble_data, scenario_grid,
    ta_col, pr_col,
    scale, area_weight, support,
    verbose, tag) {

  scale <- match.arg(scale, c("none", "global", "by_group"))
  area_weight <- match.arg(area_weight, c("regular", "none"))
  support <- .scenprob_validate_support(support)

  grid_ta <- scenario_grid[[ta_col]]
  grid_pr <- scenario_grid[[pr_col]]
  n_grid <- nrow(scenario_grid)

  dT <- .scenprob_grid_step(grid_ta)
  dP <- .scenprob_grid_step(grid_pr)

  support_mask <- .scenprob_support_mask(grid_ta, grid_pr, support)

  scale_none <- list(mean = 0, sd = 1, enabled = FALSE)

  scale_global_ta <- scale_none
  scale_global_pr <- scale_none
  grid_ta_global <- grid_ta
  grid_pr_global <- grid_pr
  dT_global <- dT
  dP_global <- dP

  if (scale == "global") {
    scale_global_ta <- .scenprob_scale_params(ensemble_data[[ta_col]], "ta", " (global)", verbose, tag)
    scale_global_pr <- .scenprob_scale_params(ensemble_data[[pr_col]], "pr", " (global)", verbose, tag)
    grid_ta_global <- .scenprob_apply_scale(grid_ta, scale_global_ta)
    grid_pr_global <- .scenprob_apply_scale(grid_pr, scale_global_pr)
    dT_global <- .scenprob_scale_step(dT, scale_global_ta)
    dP_global <- .scenprob_scale_step(dP, scale_global_pr)
  }

  area_weight_global <- .scenprob_area_weight_guard(area_weight, dT_global, dP_global, verbose, tag)

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

.scenprob_group_scale_context <- function(
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
    scale_ta <- .scenprob_scale_params(ta, "ta", paste0(" (group=", grp, ")"), verbose, tag)
    scale_pr <- .scenprob_scale_params(pr, "pr", paste0(" (group=", grp, ")"), verbose, tag)
    grid_ta_use <- .scenprob_apply_scale(grid_ctx$grid_ta, scale_ta)
    grid_pr_use <- .scenprob_apply_scale(grid_ctx$grid_pr, scale_pr)
    dT_use <- .scenprob_scale_step(grid_ctx$dT, scale_ta)
    dP_use <- .scenprob_scale_step(grid_ctx$dP, scale_pr)
    area_weight_use <- .scenprob_area_weight_guard(area_weight, dT_use, dP_use, verbose, tag)
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

.scenprob_order_weight_cols <- function(out, scenario_grid) {
  base_cols <- names(scenario_grid)
  weight_cols <- setdiff(names(out), base_cols)
  out[, c(base_cols, sort(weight_cols)), drop = FALSE]
}


# ------------------------------------------------------------------------------
# KDE-specific helpers
# ------------------------------------------------------------------------------

.scenprob_validate_bw_vec <- function(bw_vec) {
  is.numeric(bw_vec) && length(bw_vec) == 2 &&
    all(is.finite(bw_vec)) && all(bw_vec > 0)
}

.scenprob_pick_bw_auto <- function(ta, pr, dT, dP, k, alpha, bw_min, bw_max) {
  bw_nrd <- tryCatch(
    c(MASS::bandwidth.nrd(ta), MASS::bandwidth.nrd(pr)),
    error = function(e) c(NA_real_, NA_real_)
  )

  bw_grid <- c(k[1] * dT, k[2] * dP)
  if (!is.finite(bw_grid[1])) bw_grid[1] <- 0
  if (!is.finite(bw_grid[2])) bw_grid[2] <- 0

  bw_use <- pmax(bw_grid, alpha * bw_nrd, na.rm = TRUE)

  bw_use <- pmax(bw_use, bw_min)
  bw_use <- pmin(bw_use, bw_max)

  if (!.scenprob_validate_bw_vec(bw_use)) return(NULL)
  bw_use
}

#' Vectorized KDE chunk evaluation
#'
#' Computes weighted Gaussian product-kernel density for a chunk of grid points.
#' Uses matrix operations instead of explicit loops for improved performance.
#'
#' @param g_ta Numeric vector of grid temperature values (chunk).
#' @param g_pr Numeric vector of grid precipitation values (chunk).
#' @param ta_use Numeric vector of observation temperature values.
#' @param pr_use Numeric vector of observation precipitation values.
#' @param w_obs Numeric vector of observation weights (sum to 1).
#' @param inv_bw_ta Inverse bandwidth for temperature.
#' @param inv_bw_pr Inverse bandwidth for precipitation.
#' @return Numeric vector of density values for the chunk.
#' @keywords internal
.scenprob_kde_chunk_vectorized <- function(g_ta, g_pr, ta_use, pr_use,
                                           w_obs, inv_bw_ta, inv_bw_pr) {
  n_chunk <- length(g_ta)
  n_obs <- length(ta_use)

  # For moderate chunk sizes, use full matrix approach

  # outer() creates n_chunk x n_obs matrices of standardized differences
  u_mat <- outer(g_ta, ta_use, `-`) * inv_bw_ta
  v_mat <- outer(g_pr, pr_use, `-`) * inv_bw_pr

  # Gaussian kernel: dnorm(u) * dnorm(v) for each (grid, obs) pair
  # Then weighted sum across observations (matrix-vector product)
  kernel_vals <- stats::dnorm(u_mat) * stats::dnorm(v_mat)


  # Result: n_chunk vector of density contributions

  as.numeric(kernel_vals %*% w_obs)
}

#' Compute effective sample size from weights
#'
#' @param w Numeric vector of weights (need not sum to 1).
#' @return Effective sample size (Kish's formula).
#' @keywords internal
.scenprob_effective_n <- function(w) {

  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}


# ------------------------------------------------------------------------------
# MVN-specific helpers (Cholesky-based)
# ------------------------------------------------------------------------------

.scenprob_wmean <- function(x, w) {
  s <- sum(w)
  if (!is.finite(s) || s <= .Machine$double.eps) return(NA_real_)
  sum(w * x) / s
}

.scenprob_wcov2 <- function(x, y, w) {
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

#' Bivariate normal density using Cholesky decomposition
#'
#' More numerically stable than direct matrix inversion, especially for
#' near-singular covariance matrices. Uses log-space computation to prevent
#' underflow for extreme values.
#'
#' @param x Matrix (n x 2) of evaluation points.
#' @param mu Numeric vector length 2, mean.
#' @param Sigma 2x2 covariance matrix.
#' @return Numeric vector of density values, or NA if decomposition fails.
#' @keywords internal
.scenprob_dmvnorm2_chol <- function(x, mu, Sigma) {
  n <- nrow(x)


  # Cholesky decomposition: Sigma = L' L where L is upper triangular

  # More stable than computing solve(Sigma) and det(Sigma) separately

  chol_result <- tryCatch(
    chol(Sigma),
    error = function(e) NULL
  )


  if (is.null(chol_result)) {
    return(rep(NA_real_, n))
  }

  L <- chol_result  # Upper triangular: Sigma = t(L) %*% L

  # Log determinant: det(Sigma) = det(L' L) = det(L)^2

  # det(L) = prod(diag(L)), so log(det(Sigma)) = 2 * sum(log(diag(L)))
  log_det <- 2 * sum(log(diag(L)))

  # Check for degeneracy

  if (!is.finite(log_det)) {
    return(rep(NA_real_, n))
  }

  # Center the data
  xc <- cbind(x[, 1] - mu[1], x[, 2] - mu[2])

  # Solve L' z = xc' for z using backsolve

  # Then quadratic form = ||z||^2 for each observation
  # This is equivalent to xc %*% solve(Sigma) %*% t(xc) diagonal
  z <- backsolve(L, t(xc), transpose = TRUE)
  qf <- colSums(z^2)

  # Log density: -0.5 * (k*log(2*pi) + log|Sigma| + quadratic form)
  # where k = 2 for bivariate

  log_dens <- -0.5 * (2 * log(2 * pi) + log_det + qf)

  # Exponentiate, handling potential underflow gracefully
  # (exp of very negative numbers -> 0, which is fine)
  exp(log_dens)
}

#' Legacy bivariate normal density (direct inversion)
#'
#' Kept for backward compatibility and comparison. Uses direct matrix
#' inversion which can be less stable for ill-conditioned matrices.
#'
#' @param x Matrix (n x 2) of evaluation points.
#' @param mu Numeric vector length 2, mean.
#' @param Sigma 2x2 covariance matrix.
#' @return Numeric vector of density values.
#' @keywords internal
.scenprob_dmvnorm2_legacy <- function(x, mu, Sigma) {
  invS <- tryCatch(solve(Sigma), error = function(e) matrix(NA_real_, 2, 2))
  detS <- tryCatch(det(Sigma), error = function(e) NA_real_)

  if (!is.finite(detS) || detS <= .Machine$double.eps) {
    return(rep(NA_real_, nrow(x)))
  }
  if (any(!is.finite(invS))) {
    return(rep(NA_real_, nrow(x)))
  }

  xc1 <- x[, 1] - mu[1]
  xc2 <- x[, 2] - mu[2]

  qf <- invS[1, 1] * xc1 * xc1 +
    2 * invS[1, 2] * xc1 * xc2 +
    invS[2, 2] * xc2 * xc2

  norm_const <- 1 / (2 * pi * sqrt(detS))
  norm_const * exp(-0.5 * qf)
}

#' Validate covariance matrix for MVN fitting
#'
#' Checks for finite values, positive diagonal, and non-singularity.
#'
#' @param Sigma 2x2 covariance matrix.
#' @param grp Group label for error messages.
#' @return List with `valid` (logical) and `reason` (character or NULL).
#' @keywords internal
.scenprob_validate_sigma <- function(Sigma, grp = "") {
  if (any(!is.finite(Sigma))) {
    return(list(valid = FALSE, reason = "covariance matrix contains non-finite values"))
  }

  if (any(diag(Sigma) <= .Machine$double.eps)) {
    return(list(valid = FALSE, reason = "near-zero variance in fitted covariance"))
  }

  # Test Cholesky decomposition as singularity check
  chol_test <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_test)) {
    return(list(valid = FALSE, reason = "singular or non-positive-definite covariance"))
  }

  list(valid = TRUE, reason = NULL)
}


# ==============================================================================
# Main KDE function
# ==============================================================================

#' Calculate Climate-Informed Scenario Weights (KDE; auto bandwidth; optional family weighting)
#'
#' @description
#' For each group (e.g., SSP), estimates a smooth 2D kernel density surface over
#' (temperature change, precipitation change) using a Gaussian product-kernel KDE,
#' evaluates it on a user-provided stress-test grid, and returns per-group weights.
#'
#' Note that the KDE output `dens` is a continuous density, not a probability.
#' When `normalize = TRUE`, the function returns a discrete probability mass
#' function (PMF) over the provided grid points. Without area weighting, that PMF
#' depends on grid resolution. Using `area_weight = "regular"` applies regular-grid
#' cell areas before normalization, yielding results that are invariant to grid
#' resolution for regular grids (i.e., constant `dT * dP` cell area).
#'
#' If `bw = NULL`, bandwidths are chosen automatically per group using a grid-aware rule:
#' \deqn{bw = max(k * grid\_step,\ \alpha * bw\_nrd)}
#' where `grid_step` is the median spacing in `scenario_grid` for each dimension,
#' `bw_nrd` is a data-driven bandwidth estimate for each dimension, `k` controls
#' smoothing relative to grid spacing, and `alpha` scales the data-driven bandwidth.
#'
#' KDE is sensitive to variable units and relative axis magnitudes. When `scale` is
#' enabled, the KDE is computed entirely in standardized space: ensemble observations
#' and grid coordinates are centered/scaled, auto bandwidth selection happens in that
#' scaled space, and the grid-step floor uses scaled grid steps. If the user supplies
#' `bw`, it is interpreted in original units and converted to scaled units internally
#' by dividing by the corresponding standard deviation (for the chosen scaling mode).
#'
#' Optionally, observations can be downweighted by model family (near-duplicate control)
#' so that each family contributes equal total influence within each group.
#'
#' Gaussian KDE does not respect hard boundaries, so density can "leak" outside
#' plausible ranges. The `support` argument is a pragmatic truncation guard that
#' hard-zeros grid points outside user-specified bounds; it does not implement
#' boundary-corrected kernels.
#'
#' @details
#' ## Typical use cases
#' - **Stress-test weighting within scenario families**: You have an ensemble of
#'   projected changes (e.g., from multiple GCM/RCM runs) for each scenario group
#'   and you want a smooth, nonparametric weighting over a pre-defined stress-test
#'   grid of \eqn{\Delta T} and \eqn{\Delta P} combinations.
#' - **Scenario prioritization for downstream simulation**: You want a compact set
#'   of grid states to sample more frequently (high KDE mass) and deprioritize
#'   implausible combinations (low mass), while retaining group separation.
#' - **Near-duplicate control**: You suspect clusters driven by correlated model
#'   lineages; family weighting reduces over-representation of large families.
#'
#' ## Strategies for choosing key parameters
#' ### 1) `scenario_grid` design and `area_weight`
#' - If your grid is **regular** (constant spacing in both dimensions), prefer
#'   `area_weight = "regular"` with `normalize = TRUE`. This approximates an
#'   integral over the continuous KDE and makes results stable to grid refinement.
#' - If your grid is **irregular** (uneven spacing), `area_weight = "regular"`
#'   is not appropriate. This function does not implement irregular cell areas,
#'   so you should either:
#'   (i) use `area_weight = "none"` (accept grid-dependent PMF), or
#'   (ii) convert to a regular grid upstream, or
#'   (iii) implement external area weights and post-multiply before normalization.
#'
#' Practical rule:
#' - When comparing weights across alternative grids or resolutions, do not use
#'   `area_weight = "none"` unless you intentionally want grid-dependent masses.
#'
#' ### 2) Bandwidth selection: `bw` vs (`k`, `alpha`, `bw_min`, `bw_max`)
#' Bandwidth governs the bias–variance trade-off:
#' - Too small: weights become spiky (overfit), sensitive to sampling noise and
#'   duplicate runs.
#' - Too large: weights become overly diffuse, blurring distinctions between
#'   plausible and implausible regions.
#'
#' Recommended workflows:
#' - **Default exploratory run**: `bw = NULL`, keep defaults (`k`, `alpha`) and
#'   inspect resulting weight surfaces. Use this when you do not have strong prior
#'   beliefs about smoothing scale.
#' - **Grid-limited smoothing**: Increase `k` when your grid spacing is coarse and
#'   you want at least a few grid steps of smoothing (avoid single-cell peaks).
#'   Decrease `k` when the grid is dense and you want more local structure.
#' - **Data-limited smoothing**: Increase `alpha` when sample size per group is
#'   small or noisy (more smoothing). Decrease `alpha` when groups are large and
#'   you want to preserve multimodality.
#'
#' Use bounds to stabilize edge cases:
#' - `bw_min` prevents degenerately small bandwidths (spikes) in groups with
#'   near-collinearity or low variance.
#' - `bw_max` caps oversmoothing when `bandwidth.nrd()` becomes large due to
#'   broad dispersion or outliers.
#'
#' Recommended starting points (then adjust based on diagnostics):
#' - Coarse grid or discrete stress-test states: raise `k` (e.g., 2–4).
#' - Very dense grid: lower `k` (e.g., 1–1.5) to avoid excessive diffusion.
#' - Small n per group (< ~20): raise `alpha` (e.g., 1.5–2).
#' - Large n per group (> ~100): `alpha` near 1 is usually sufficient.
#'
#' If you need strict comparability across groups, supply a fixed `bw` (in original
#' units) so all groups use identical smoothing. This is useful when differences
#' in auto-selected bandwidths are confounded with differences in group sample sizes.
#'
#' ### 3) Scaling: `scale = "none" | "global" | "by_group"`
#' KDE geometry depends on relative axis scaling. Choose scaling based on your
#' objective:
#' - `"none"`: Use when both variables are already commensurate or intentionally
#'   expressed in comparable units (e.g., both standardized anomalies), or when
#'   you want bandwidths interpretable directly in physical units without any
#'   implicit rescaling.
#' - `"global"`: Use when units/magnitudes differ (common for \eqn{\Delta T} in °C
#'   vs \eqn{\Delta P} in percent) and you want a consistent kernel geometry across
#'   groups. This is the most defensible default for multi-group comparison.
#' - `"by_group"`: Use when within-group dispersion differs substantially and you
#'   want KDE to adapt to each group's spread in standardized space. This increases
#'   within-group resolution but reduces comparability of absolute smoothing scales
#'   across groups (because each group uses its own SDs).
#'
#' Practical cautions:
#' - If you supply `bw`, it is always interpreted in original units. Under scaling,
#'   it will be internally converted using SDs; therefore, group-to-group effective
#'   smoothing differs under `"by_group"` even with the same `bw` in original units.
#'
#' ### 4) Family downweighting: `use_family_weights`, `family_col`
#' Use family weighting when you have multiple realizations that are not
#' statistically independent (e.g., many closely related model variants):
#' - With `use_family_weights = TRUE`, each family contributes equal total weight
#'   within a group, and members within a family split that family share evenly.
#' - Use this when you want to treat families as the effective sample units.
#'
#' Do not use family weighting if:
#' - Your "family" labels are noisy/ambiguous, or
#' - You explicitly want the ensemble to reflect the available run counts.
#'
#' ### 5) Truncation guard: `support`
#' Use `support` to prevent allocating weight to physically or programmatically
#' irrelevant parts of the grid (e.g., negative precipitation-change bounds if your
#' stress test excludes drying). This is not boundary correction; it is a hard mask.
#'
#' Strategies:
#' - Set `support` to the plausible envelope implied by your scenario definition
#'   or by expert judgment (e.g., exclude combinations outside a policy-relevant
#'   design space).
#' - Prefer masking on the *grid* (as done here) rather than trimming the ensemble,
#'   unless you want to change the KDE fit itself.
#'
#' ### 6) Data sufficiency and skipping: `min_samples`
#' Groups with too few complete observations or near-zero variance in either
#' dimension are skipped. Set `min_samples` based on how stable you need KDE
#' estimates to be:
#' - Increase `min_samples` when you require robust multimodality and stable tails.
#' - Decrease `min_samples` only if you accept noisy, potentially spurious surfaces.
#'
#' Operationally, examine `attr(out, "skipped_groups")` and decide whether to
#' (i) merge sparse groups, (ii) expand the ensemble, or (iii) treat skipped groups
#' as requiring manual weights.
#'
#' ### 7) Performance and memory: `chunk_size`
#' The KDE evaluation scales as O(n_obs * n_grid) per group. `chunk_size` controls
#' the size of grid blocks used during evaluation:
#' - Larger `chunk_size` reduces overhead but increases peak memory.
#' - Smaller `chunk_size` reduces peak memory and may be needed for large grids,
#'   at the cost of more looping overhead.
#'
#' Practical guidance:
#' - For grids in the 10^4–10^5 range, start at 5,000–20,000.
#' - If you see memory pressure, reduce `chunk_size` before changing the grid.
#'
#' ### 8) Diagnostics: `diagnostics`
#' When `diagnostics = TRUE`, additional attributes are attached to the output:
#' - `bandwidth_used`: Named list of bandwidths (in working space) per group.
#' - `effective_sample_size`: Named list of effective n per group (after family weighting).
#' - `scaling_params`: Scaling parameters used (global or per-group).
#'
#' @param ensemble_data
#'   Data frame containing the ensemble of scenario change pairs used to fit the
#'   KDE within each group. Must contain `ta_col`, `pr_col`, and `group_col`.
#'   Rows with non-finite values in either variable are dropped within each group.
#'
#' @param scenario_grid
#'   Data frame defining the stress-test design space where KDE is evaluated.
#'   Must contain columns named `ta_col` and `pr_col`. Each row is one grid point.
#'   The output retains this grid and appends one weight column per group.
#'
#' @param pr_col
#'   Character scalar giving the precipitation-change column name in both
#'   `ensemble_data` and `scenario_grid`. Use a consistent representation
#'   (e.g., percent change, absolute change) across both inputs.
#'
#' @param ta_col
#'   Character scalar giving the temperature-change column name in both
#'   `ensemble_data` and `scenario_grid` (e.g., °C change).
#'
#' @param group_col
#'   Character scalar naming the grouping variable in `ensemble_data` that defines
#'   separate KDE fits (e.g., `"scenario"`, `"ssp"`, `"pathway"`). Each unique value
#'   becomes one output weight column (made syntactically valid via `make.names()`).
#'
#' @param bw
#'   Optional numeric length-2 vector giving fixed bandwidths in **original units**:
#'   `c(bw_ta, bw_pr)`. If supplied, the same bandwidth is used for all groups.
#'   If `NULL` (default), bandwidths are auto-selected per group using `k`, `alpha`,
#'   and the grid step. Under `scale != "none"`, fixed `bw` is internally converted
#'   to scaled units using the relevant SDs.
#'
#' @param k
#'   Numeric length-2 vector of **grid-step multipliers** used only when `bw = NULL`.
#'   The product `k * grid_step` defines a floor on the bandwidth in each dimension
#'   to prevent bandwidths that are too small relative to grid resolution. Increase
#'   `k` to enforce smoother, less spiky weights on coarse grids; decrease `k` to
#'   allow finer structure on dense grids.
#'
#' @param alpha
#'   Numeric scalar multiplier applied to the data-driven bandwidth estimates
#'   (`bandwidth.nrd`) when `bw = NULL`. Values > 1 increase smoothing (useful for
#'   small or noisy groups); values close to 1 keep the default data-driven scale.
#'
#' @param bw_min
#'   Numeric length-2 vector giving minimum allowable bandwidths (in original units
#'   if `scale = "none"`, otherwise converted to scaled units internally). Use this
#'   to prevent near-zero bandwidths that yield unstable, highly localized weights.
#'
#' @param bw_max
#'   Numeric length-2 vector giving maximum allowable bandwidths (in original units
#'   if `scale = "none"`, otherwise converted to scaled units internally). Use this
#'   to cap oversmoothing (e.g., in the presence of outliers or broad dispersion).
#'
#' @param min_samples
#'   Integer scalar giving the minimum number of complete (finite) observations
#'   required per group to fit KDE. Groups failing this threshold are skipped and
#'   listed in `attr(out, "skipped_groups")`. Must be >= 3; recommended higher if
#'   you need stable multimodal shapes.
#'
#' @param normalize
#'   Logical scalar. If `TRUE` (default), normalize weights within each group so the
#'   resulting grid PMF sums to 1 (after optional `area_weight` and `support` masking).
#'   If `FALSE`, returns raw KDE densities evaluated on the grid (not comparable across
#'   grids without further processing).
#'
#' @param area_weight
#'   Character scalar controlling whether grid-cell area is applied before normalization.
#'   - `"regular"`: multiplies densities by `dT * dP` (median grid step products) prior
#'     to normalization, approximating a continuous integral on a **regular** grid and
#'     making normalized weights stable to grid refinement.
#'   - `"none"`: normalizes raw pointwise densities; the resulting PMF depends on the
#'     grid resolution and spacing.
#'   If grid steps are invalid or degenerate, `"regular"` is automatically disabled.
#'
#' @param scale
#'   Character scalar controlling whether KDE is computed in standardized space:
#'   - `"none"`: no scaling; KDE geometry reflects original units.
#'   - `"global"`: scale using overall mean/SD across all groups (consistent geometry).
#'   - `"by_group"`: scale within each group (adaptive geometry, weaker cross-group
#'     comparability). Scaling is disabled for a dimension if its SD is non-finite or
#'     near zero in the relevant scope.
#'
#' @param use_family_weights
#'   Logical scalar. If `TRUE`, apply observation weights so each `family_col` level
#'   contributes equal total influence within each group (near-duplicate control).
#'   If `FALSE` (default), all observations contribute equally within each group.
#'
#' @param family_col
#'   Character scalar naming the column in `ensemble_data` that identifies model family
#'   membership for `use_family_weights = TRUE`. Missing values are treated as a single
#'   `"NA_FAMILY"` group for weighting purposes.
#'
#' @param support
#'   Optional named list specifying hard bounds in original units used to mask the
#'   **grid** (not the input ensemble). Provide `ta = c(min, max)` and/or
#'   `pr = c(min, max)`. Grid points outside the specified bounds receive zero weight
#'   before normalization. This is a truncation guard, not boundary-corrected KDE.
#'
#' @param chunk_size
#'   Integer scalar controlling the number of grid points processed per chunk during
#'   KDE evaluation. Used to limit peak memory for large grids. Performance scales
#'   approximately with O(n_obs * n_grid) per group; `chunk_size` trades off overhead
#'   versus peak memory.
#'
#' @param diagnostics
#'   Logical scalar. If `TRUE`, attach additional attributes with per-group diagnostic
#'   information: `bandwidth_used`, `effective_sample_size`, and `scaling_params`.
#'   Default is `FALSE`.
#'
#' @param verbose
#'   Logical scalar. If `TRUE`, prints per-group messages about bandwidth selection,
#'   scaling disablement, and skipped groups. Use `FALSE` for silent batch workflows.
#'
#' @return
#' Data frame identical to `scenario_grid` plus one weight column per group
#' (alphabetically ordered by group name after `make.names()`).
#'
#' Attributes:
#' - `"skipped_groups"`: Character vector of groups skipped due to insufficient data
#'   or degeneracy.
#'
#' When `diagnostics = TRUE`, additional attributes:
#' - `"bandwidth_used"`: Named list with bandwidth vectors (in working space) per group.
#' - `"effective_sample_size"`: Named list with effective n per group.
#' - `"scaling_params"`: List with global or per-group scaling parameters.
#'
#' @examples
#' # --- Example 1: Basic usage with auto bandwidth (two groups) ---
#' set.seed(1)
#' ensemble_data <- data.frame(
#'   scenario = rep(c("SSP1", "SSP2"), each = 30),
#'   tavg     = c(rnorm(30, mean = 1.5, sd = 0.4), rnorm(30, mean = 2.5, sd = 0.5)),
#'   prcp     = c(rnorm(30, mean = 5,   sd = 1.0), rnorm(30, mean = 0,   sd = 1.2))
#' )
#'
#' scenario_grid <- expand.grid(
#'   tavg = seq(0.5, 3.5, by = 0.25),
#'   prcp = seq(-3,  8,   by = 0.5)
#' )
#'
#' w <- estimate_scenario_probs_kde(
#'   ensemble_data  = ensemble_data,
#'   scenario_grid  = scenario_grid,
#'   group_col      = "scenario",
#'   ta_col         = "tavg",
#'   pr_col         = "prcp",
#'   bw             = NULL,
#'   normalize      = TRUE,
#'   verbose        = FALSE
#' )
#'
#' # Each group column sums to 1 over the grid (when normalize = TRUE)
#' colSums(w[, c("SSP1", "SSP2")])
#'
#' # --- Example 2: With diagnostics output ---
#' w_diag <- estimate_scenario_probs_kde(
#'   ensemble_data  = ensemble_data,
#'   scenario_grid  = scenario_grid,
#'   diagnostics    = TRUE,
#'   verbose        = FALSE
#' )
#' attr(w_diag, "bandwidth_used")
#' attr(w_diag, "effective_sample_size")
#'
#' @importFrom MASS bandwidth.nrd
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
estimate_scenario_probs_kde <- function(
    ensemble_data,
    scenario_grid,
    pr_col = "prcp",
    ta_col = "tavg",
    group_col = "scenario",
    bw = NULL,
    k = c(1.5, 2.0),
    alpha = 1.0,
    bw_min = c(0, 0),
    bw_max = c(Inf, Inf),
    min_samples = 5L,
    normalize = TRUE,
    area_weight = c("regular", "none"),
    scale = c("none", "global", "by_group"),
    use_family_weights = FALSE,
    family_col = "model_family",
    support = NULL,
    chunk_size = 5000L,
    diagnostics = FALSE,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation (keep strict)
  # ---------------------------------------------------------------------------

  required_cols <- c(ta_col, pr_col, group_col)
  missing_cols <- setdiff(required_cols, names(ensemble_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ensemble_data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!all(c(ta_col, pr_col) %in% names(scenario_grid))) {
    stop("Both temperature and precipitation columns must exist in scenario_grid.",
         call. = FALSE)
  }

  if (nrow(scenario_grid) == 0) stop("scenario_grid cannot be empty.", call. = FALSE)

  if (!is.logical(normalize) || length(normalize) != 1) {
    stop("normalize must be a single logical.", call. = FALSE)
  }

  if (!is.logical(diagnostics) || length(diagnostics) != 1) {
    stop("diagnostics must be a single logical.", call. = FALSE)
  }


  area_weight <- match.arg(area_weight)
  scale <- match.arg(scale)

  if (!is.numeric(min_samples) || length(min_samples) != 1 || min_samples < 3) {
    stop("min_samples must be a single integer >= 3.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  if (!is.null(bw) && !.scenprob_validate_bw_vec(bw)) {
    stop("If provided, bw must be numeric length-2 with finite values > 0: c(bw_ta, bw_pr).",
         call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 2 || any(!is.finite(k)) || any(k <= 0)) {
    stop("k must be numeric length-2 with finite values > 0.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a single finite numeric value > 0.", call. = FALSE)
  }

  if (!is.numeric(bw_min) || length(bw_min) != 2 || any(!is.finite(bw_min)) || any(bw_min < 0)) {
    stop("bw_min must be numeric length-2 with finite values >= 0.", call. = FALSE)
  }

  if (!is.numeric(bw_max) || length(bw_max) != 2 || any(is.na(bw_max)) || any(is.nan(bw_max))) {
    stop("bw_max must be numeric length-2 with non-missing values.", call. = FALSE)
  }
  if (any(bw_max <= 0)) stop("bw_max must be > 0 (use Inf for no upper bound).", call. = FALSE)
  if (any(bw_max <= bw_min)) stop("bw_max must be greater than bw_min for both dimensions.", call. = FALSE)

  if (use_family_weights && !(family_col %in% names(ensemble_data))) {
    stop("use_family_weights=TRUE but family_col not found in ensemble_data: ", family_col,
         call. = FALSE)
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 100) {
    stop("chunk_size must be a single integer >= 100.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  # ---------------------------------------------------------------------------
  # Grid context (shared)
  # ---------------------------------------------------------------------------

  grid_ctx <- .scenprob_prepare_grid_context(
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
  # Prepare output, diagnostics storage, and loop
  # ---------------------------------------------------------------------------

  out <- scenario_grid
  groups <- sort(unique(as.character(ensemble_data[[group_col]])))
  skipped <- character(0)

  # Diagnostics storage
  diag_bandwidth <- list()
  diag_eff_n <- list()
  diag_scaling <- list()

  for (grp in groups) {

    df <- dplyr::filter(ensemble_data, .data[[group_col]] == grp)

    ta <- df[[ta_col]]
    pr <- df[[pr_col]]

    ok <- is.finite(ta) & is.finite(pr)
    ta <- ta[ok]
    pr <- pr[ok]

    if (length(ta) < min_samples) {
      warning("Skipping group '", grp, "' (need at least ", min_samples, " complete observations).",
              call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    if (stats::sd(ta) < .Machine$double.eps || stats::sd(pr) < .Machine$double.eps) {
      warning("Skipping group '", grp, "' (near-zero variance in temperature or precipitation).",
              call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    # Group scaling context (shared)
    sc_ctx <- .scenprob_group_scale_context(
      ta = ta, pr = pr,
      grid_ctx = grid_ctx,
      grp = grp,
      verbose = verbose,
      tag = "KDE",
      area_weight = area_weight
    )

    ta_use <- .scenprob_apply_scale(ta, sc_ctx$scale_ta)
    pr_use <- .scenprob_apply_scale(pr, sc_ctx$scale_pr)

    bw_min_use <- .scenprob_scale_bw(bw_min, sc_ctx$scale_ta, sc_ctx$scale_pr)
    bw_max_use <- .scenprob_scale_bw(bw_max, sc_ctx$scale_ta, sc_ctx$scale_pr)

    # Bandwidth choice
    bw_use <- bw
    if (is.null(bw_use)) {
      bw_use <- .scenprob_pick_bw_auto(
        ta = ta_use, pr = pr_use,
        dT = sc_ctx$dT_use, dP = sc_ctx$dP_use,
        k = k, alpha = alpha,
        bw_min = bw_min_use, bw_max = bw_max_use
      )

      if (is.null(bw_use)) {
        warning("Skipping group '", grp, "' (auto bandwidth selection failed).", call. = FALSE)
        skipped <- c(skipped, grp)
        next
      }
    } else {
      bw_use <- .scenprob_scale_bw(bw_use, sc_ctx$scale_ta, sc_ctx$scale_pr)
    }

    if (isTRUE(verbose)) {
      msg <- if (is.null(bw)) "auto" else "user"
      bw_label <- if (scale == "none") "bw" else "bw_scaled"
      message(
        sprintf(
          "[KDE] Group=%s | %s_%s = (ta=%.3f, pr=%.3f)",
          grp, bw_label, msg, bw_use[1], bw_use[2]
        )
      )
    }

    bw_ta <- bw_use[1]
    bw_pr <- bw_use[2]

    # Observation weights (uniform or family-weighted)
    if (use_family_weights) {
      fam <- df[[family_col]][ok]
      w_obs <- .scenprob_compute_family_weights(fam)
      w_obs <- w_obs / sum(w_obs)
    } else {
      w_obs <- rep(1 / length(ta), length(ta))
    }

    # Store diagnostics
    if (diagnostics) {
      diag_bandwidth[[grp]] <- c(ta = bw_ta, pr = bw_pr)
      diag_eff_n[[grp]] <- .scenprob_effective_n(w_obs * length(w_obs))  # unnormalized for ESS
      diag_scaling[[grp]] <- list(
        ta = sc_ctx$scale_ta,
        pr = sc_ctx$scale_pr
      )
    }

    # -------------------------------------------------------------------------
    # Vectorized KDE evaluation
    # -------------------------------------------------------------------------
    n_obs <- length(ta_use)
    dens <- numeric(n_grid)

    inv_bw_ta <- 1 / bw_ta
    inv_bw_pr <- 1 / bw_pr
    norm_const <- inv_bw_ta * inv_bw_pr

    idx_starts <- seq.int(1L, n_grid, by = chunk_size)

    for (s in idx_starts) {
      e <- as.integer(min(s + chunk_size - 1L, n_grid))

      g_ta <- sc_ctx$grid_ta_use[s:e]
      g_pr <- sc_ctx$grid_pr_use[s:e]

      # Vectorized chunk computation
      chunk_dens <- .scenprob_kde_chunk_vectorized(
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
      weights <- dens * (sc_ctx$dT_use * sc_ctx$dP_use)
    }

    if (!is.null(support_mask)) weights[!support_mask] <- 0
    if (normalize) weights <- .scenprob_normalize_vec(weights)

    out[[make.names(grp)]] <- weights
  }

  out <- .scenprob_order_weight_cols(out, scenario_grid)

  # ---------------------------------------------------------------------------
  # Attach attributes
  # ---------------------------------------------------------------------------

  if (length(skipped) > 0) {
    if (isTRUE(verbose)) message("Skipped ", length(skipped), " group(s): ", paste(skipped, collapse = ", "))
    attr(out, "skipped_groups") <- skipped
  } else {
    attr(out, "skipped_groups") <- character(0)
  }

  if (diagnostics) {
    attr(out, "bandwidth_used") <- diag_bandwidth
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
# Main MVN function
# ==============================================================================

#' Calculate Climate-Informed Scenario Weights (MVN; optional robust fit; optional family weighting)
#'
#' @description
#' For each group (e.g., SSP), fits a bivariate normal distribution to
#' (temperature change, precipitation change) and evaluates the fitted density
#' on a user-provided stress-test grid. Returns per-group weights.
#'
#' As with KDE, the MVN density is a continuous density, not a probability. When
#' `normalize = TRUE`, the function returns a discrete probability mass function
#' (PMF) over the provided grid points. Without area weighting, that PMF depends
#' on grid resolution. Using `area_weight = "regular"` applies regular-grid cell
#' areas before normalization, yielding results that are invariant to grid
#' resolution for regular grids (i.e., constant `dT * dP` cell area).
#'
#' MVN is a parametric analogue to `estimate_scenario_probs_kde()`.
#' It is most appropriate when the ensemble cloud is reasonably elliptical
#' and unimodal within each group.
#'
#' This implementation uses Cholesky decomposition for numerical stability when
#' evaluating the bivariate normal density, which is more robust than direct
#' matrix inversion for near-singular covariance matrices.
#'
#' @param ensemble_data Data frame with ensemble points used to fit per-group MVN.
#' @param scenario_grid Data frame with grid points where weights are evaluated.
#' @param pr_col Character scalar. Precipitation change column name.
#' @param ta_col Character scalar. Temperature change column name.
#' @param group_col Character scalar. Grouping column in `ensemble_data`.
#' @param robust Logical scalar. If `TRUE`, use robust mean/cov via `MASS::cov.trob()`.
#' @param min_samples Integer scalar. Minimum complete observations required per group.
#'   Recommended: >= 5 for robust fits.
#' @param normalize Logical scalar. If `TRUE`, normalize weights within each group.
#' @param area_weight Character scalar: `"regular"` or `"none"`.
#' @param scale Character scalar: `"none"`, `"global"`, or `"by_group"`.
#' @param use_family_weights Logical scalar. If `TRUE`, each family contributes
#'   equal total influence within each group when estimating mean/cov.
#' @param family_col Character scalar naming the family column (if used).
#' @param support Optional named list with `ta` and/or `pr` ranges to hard-mask grid.
#' @param diagnostics Logical scalar. If `TRUE`, attach additional attributes with
#'   per-group diagnostic information: `mvn_params`, `effective_sample_size`, and
#'   `scaling_params`. Default is `FALSE`.
#' @param verbose Logical scalar. If `TRUE`, prints per-group messages.
#'
#' @return
#' Data frame identical to `scenario_grid` plus one weight column per group
#' (alphabetically ordered by group name after `make.names()`).
#'
#' Attributes:
#' - `"skipped_groups"`: Character vector of groups skipped due to insufficient data
#'   or degeneracy.
#'
#' When `diagnostics = TRUE`, additional attributes:
#' - `"mvn_params"`: Named list with `mu` and `Sigma` (in working space) per group.
#' - `"effective_sample_size"`: Named list with effective n per group.
#' - `"scaling_params"`: List with global or per-group scaling parameters.
#' - `"fit_method"`: Character indicating "robust" or "classical" per group.
#'
#' @examples
#' set.seed(42)
#' ensemble_data <- data.frame(
#'   scenario = rep(c("SSP1", "SSP2"), each = 50),
#'   tavg     = c(rnorm(50, 1.5, 0.4), rnorm(50, 2.5, 0.5)),
#'   prcp     = c(rnorm(50, 5, 1.0), rnorm(50, 0, 1.2))
#' )
#'
#' scenario_grid <- expand.grid(
#'   tavg = seq(0, 4, by = 0.2),
#'   prcp = seq(-4, 10, by = 0.5)
#' )
#'
#' # Basic usage
#' w <- estimate_scenario_probs_mvnorm(
#'   ensemble_data = ensemble_data,
#'   scenario_grid = scenario_grid,
#'   verbose = FALSE
#' )
#' colSums(w[, c("SSP1", "SSP2")])
#'
#' # With diagnostics
#' w_diag <- estimate_scenario_probs_mvnorm(
#'   ensemble_data = ensemble_data,
#'   scenario_grid = scenario_grid,
#'   diagnostics = TRUE,
#'   verbose = FALSE
#' )
#' attr(w_diag, "mvn_params")
#'
#' @importFrom MASS cov.trob
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
estimate_scenario_probs_mvnorm <- function(
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
    use_family_weights = FALSE,
    family_col = "model_family",
    support = NULL,
    diagnostics = FALSE,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  required_cols <- c(ta_col, pr_col, group_col)
  missing_cols <- setdiff(required_cols, names(ensemble_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ensemble_data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!all(c(ta_col, pr_col) %in% names(scenario_grid))) {
    stop("Both temperature and precipitation columns must exist in scenario_grid.",
         call. = FALSE)
  }

  if (nrow(scenario_grid) == 0) stop("scenario_grid cannot be empty.", call. = FALSE)

  if (!is.logical(robust) || length(robust) != 1) stop("robust must be a single logical.", call. = FALSE)
  if (!is.logical(normalize) || length(normalize) != 1) stop("normalize must be a single logical.", call. = FALSE)
  if (!is.logical(diagnostics) || length(diagnostics) != 1) stop("diagnostics must be a single logical.", call. = FALSE)

  if (!is.numeric(min_samples) || length(min_samples) != 1 || min_samples < 3) {
    stop("min_samples must be a single integer >= 3.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  area_weight <- match.arg(area_weight)
  scale <- match.arg(scale)

  if (use_family_weights && !(family_col %in% names(ensemble_data))) {
    stop("use_family_weights=TRUE but family_col not found in ensemble_data: ", family_col,
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Grid context (shared)
  # ---------------------------------------------------------------------------

  grid_ctx <- .scenprob_prepare_grid_context(
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
  # Prepare output, diagnostics storage, and loop
  # ---------------------------------------------------------------------------

  out <- scenario_grid
  groups <- sort(unique(as.character(ensemble_data[[group_col]])))
  skipped <- character(0)

  # Diagnostics storage
  diag_mvn_params <- list()
  diag_eff_n <- list()
  diag_scaling <- list()
  diag_fit_method <- list()

  for (grp in groups) {

    df <- dplyr::filter(ensemble_data, .data[[group_col]] == grp)

    ta <- df[[ta_col]]
    pr <- df[[pr_col]]

    ok <- is.finite(ta) & is.finite(pr)
    ta <- ta[ok]
    pr <- pr[ok]

    if (length(ta) < min_samples) {
      warning("Skipping group '", grp, "' (need at least ", min_samples, " complete observations).",
              call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    if (stats::sd(ta) < .Machine$double.eps || stats::sd(pr) < .Machine$double.eps) {
      warning("Skipping group '", grp, "' (near-zero variance in temperature or precipitation).",
              call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    sc_ctx <- .scenprob_group_scale_context(
      ta = ta, pr = pr,
      grid_ctx = grid_ctx,
      grp = grp,
      verbose = verbose,
      tag = "MVN",
      area_weight = area_weight
    )

    ta_use <- .scenprob_apply_scale(ta, sc_ctx$scale_ta)
    pr_use <- .scenprob_apply_scale(pr, sc_ctx$scale_pr)

    # Estimation weights (only used for classical fit)
    if (use_family_weights) {
      fam <- df[[family_col]][ok]
      w_obs <- .scenprob_compute_family_weights(fam)
      if (is.finite(sum(w_obs)) && sum(w_obs) > .Machine$double.eps) {
        w_obs <- w_obs / sum(w_obs)
      } else {
        w_obs <- rep(1 / length(ta_use), length(ta_use))
      }
    } else {
      w_obs <- rep(1 / length(ta_use), length(ta_use))
    }

    # Store effective sample size for diagnostics
    if (diagnostics) {
      diag_eff_n[[grp]] <- .scenprob_effective_n(w_obs * length(w_obs))
      diag_scaling[[grp]] <- list(
        ta = sc_ctx$scale_ta,
        pr = sc_ctx$scale_pr
      )
    }

    fit <- tryCatch({
      if (robust) {
        # cov.trob does not accept external weights; robust fit ignores w_obs by design.
        m <- MASS::cov.trob(cbind(ta_use, pr_use))
        list(mu = as.numeric(m$center), Sigma = m$cov, method = "robust")
      } else {
        mu <- c(.scenprob_wmean(ta_use, w_obs), .scenprob_wmean(pr_use, w_obs))
        Sigma <- .scenprob_wcov2(ta_use, pr_use, w_obs)
        list(mu = mu, Sigma = Sigma, method = "classical")
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

    # Validate covariance matrix
    sigma_check <- .scenprob_validate_sigma(Sigma, grp)
    if (!sigma_check$valid) {
      warning("Skipping group '", grp, "' (", sigma_check$reason, ").", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    if (isTRUE(verbose)) {
      message(sprintf("[MVN] Group=%s | fit=%s | mu=(ta=%.3f, pr=%.3f)",
                      grp, fit$method, mu[1], mu[2]))
    }

    # Store diagnostics
    if (diagnostics) {
      diag_mvn_params[[grp]] <- list(mu = mu, Sigma = Sigma)
      diag_fit_method[[grp]] <- fit$method
    }

    # -------------------------------------------------------------------------
    # Cholesky-based density evaluation
    # -------------------------------------------------------------------------
    x_grid <- cbind(sc_ctx$grid_ta_use, sc_ctx$grid_pr_use)
    dens <- .scenprob_dmvnorm2_chol(x_grid, mu = mu, Sigma = Sigma)

    if (all(!is.finite(dens))) {
      warning("Skipping group '", grp, "' (density evaluation failed).", call. = FALSE)
      skipped <- c(skipped, grp)
      next
    }

    weights <- dens

    if (normalize && sc_ctx$area_weight_use == "regular") {
      weights <- weights * (sc_ctx$dT_use * sc_ctx$dP_use)
    }

    if (!is.null(support_mask)) weights[!support_mask] <- 0
    if (normalize) weights <- .scenprob_normalize_vec(weights)

    out[[make.names(grp)]] <- weights
  }

  out <- .scenprob_order_weight_cols(out, scenario_grid)

  # ---------------------------------------------------------------------------
  # Attach attributes
  # ---------------------------------------------------------------------------

  if (length(skipped) > 0) {
    if (isTRUE(verbose)) message("Skipped ", length(skipped), " group(s): ", paste(skipped, collapse = ", "))
    attr(out, "skipped_groups") <- skipped
  } else {
    attr(out, "skipped_groups") <- character(0)
  }

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
