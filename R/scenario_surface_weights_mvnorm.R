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
