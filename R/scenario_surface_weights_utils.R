# ==============================================================================
# Climate Scenario Weight Estimation
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

#' Generate unique column name (no transformation)
#' @keywords internal
.scenwgt_unique_colname <- function(grp, existing_names) {

  if (length(grp) != 1L) {
    stop("Group name must be length-1, got length ", length(grp), ".", call. = FALSE)
  }

  col_name <- as.character(grp)

  if (!nzchar(col_name)) {
    stop("Group name must be a non-empty string.", call. = FALSE)
  }

  if (col_name %in% existing_names) {
    stop(
      "Duplicate group/column name '", col_name, "'. ",
      "Group names must be unique and must not collide with existing output columns.",
      call. = FALSE
    )
  }

  col_name
}

