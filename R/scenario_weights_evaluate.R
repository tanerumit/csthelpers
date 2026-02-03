# ==============================================================================
# Surface-weights goodness-of-fit evaluator (method-agnostic)
# ==============================================================================

#' Evaluate scenario surface weights against ensemble points (log score on surface)
#'
#' @description
#' Evaluates a *discrete* surface-weight representation (scenario_surface + per-group
#' weight columns) by assigning each observation to the surface and computing the
#' log score: mean(log(p_assigned)).
#'
#' This function is method-agnostic: KDE/MVN/copula outputs can all be evaluated
#' as long as they return normalized surface weights (sum to 1) per group.
#'
#' @param surface_weights data.frame returned by compute_scenario_surface_weights_*().
#'   Must contain ta_col and pr_col plus one weight column per group.
#' @param ensemble_data data.frame with ta_col, pr_col, group_col and optional weights_col.
#' @param ta_col,pr_col,group_col column names.
#' @param weights_col optional per-row weights for scoring aggregation (nonnegative).
#' @param mapping mapping rule from point -> surface mass:
#'   "nearest" assigns probability mass of nearest surface point;
#'   "knn" uses k nearest points and a Gaussian distance kernel.
#' @param k number of neighbors for mapping="knn" (ignored for nearest).
#' @param kernel_bandwidth bandwidth (in scaled distance units) for mapping="knn".
#'   If NULL, uses median distance to k-th neighbor across points (per group).
#' @param eps probability floor to avoid -Inf when mass is zero.
#' @param scale "none" or "zscore": if "zscore", distances are computed in standardized
#'   space using ensemble_data means/sds (global).
#' @param chunk_size chunk size for distance computations.
#' @return data.frame with per-group metrics + OVERALL row.
#' @export
evaluate_scenario_surface_weights <- function(
    surface_weights,
    ensemble_data,
    ta_col = "tavg",
    pr_col = "prcp",
    group_col = "scenario",
    weights_col = NULL,
    mapping = c("nearest", "knn"),
    k = 5L,
    kernel_bandwidth = NULL,
    eps = 1e-15,
    scale = c("none", "zscore"),
    chunk_size = 2000L
) {
  mapping <- match.arg(mapping)
  scale <- match.arg(scale)

  # ---- validation ----
  if (!is.data.frame(surface_weights)) stop("surface_weights must be a data.frame.", call. = FALSE)
  if (!is.data.frame(ensemble_data)) stop("ensemble_data must be a data.frame.", call. = FALSE)

  if (!all(c(ta_col, pr_col) %in% names(surface_weights))) {
    stop("surface_weights must contain ta_col and pr_col columns.", call. = FALSE)
  }
  if (!all(c(ta_col, pr_col, group_col) %in% names(ensemble_data))) {
    stop("ensemble_data must contain ta_col, pr_col, and group_col columns.", call. = FALSE)
  }

  if (!is.null(weights_col)) {
    if (!is.character(weights_col) || length(weights_col) != 1L || !nzchar(weights_col)) {
      stop("weights_col must be NULL or a single non-empty character column name.", call. = FALSE)
    }
    if (!(weights_col %in% names(ensemble_data))) stop("weights_col not found in ensemble_data.", call. = FALSE)
    w <- ensemble_data[[weights_col]]
    if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0)) stop("weights_col must be finite and >= 0.", call. = FALSE)
    if (all(w == 0)) stop("weights_col is all zero.", call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 1) {
    stop("eps must be a single numeric in (0,1).", call. = FALSE)
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || chunk_size < 100L) {
    stop("chunk_size must be a single integer >= 100.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  if (mapping == "knn") {
    if (!is.numeric(k) || length(k) != 1L || k < 2L) stop("k must be integer >= 2 for mapping='knn'.", call. = FALSE)
    k <- as.integer(k)
  }

  # ---- identify group weight columns ----
  base_cols <- c(ta_col, pr_col)
  weight_cols <- setdiff(names(surface_weights), base_cols)

  if (length(weight_cols) == 0L) stop("surface_weights has no group weight columns.", call. = FALSE)

  # ---- common scaling for distance computations ----
  ens_ta <- ensemble_data[[ta_col]]
  ens_pr <- ensemble_data[[pr_col]]
  ok_ens <- is.finite(ens_ta) & is.finite(ens_pr) & !is.na(ensemble_data[[group_col]])
  ensemble_ok <- ensemble_data[ok_ens, , drop = FALSE]

  if (nrow(ensemble_ok) == 0L) stop("No finite ensemble points to score.", call. = FALSE)

  if (scale == "zscore") {
    mu_ta <- mean(ensemble_ok[[ta_col]])
    mu_pr <- mean(ensemble_ok[[pr_col]])
    sd_ta <- stats::sd(ensemble_ok[[ta_col]])
    sd_pr <- stats::sd(ensemble_ok[[pr_col]])
    if (!is.finite(sd_ta) || sd_ta <= .Machine$double.eps) sd_ta <- 1
    if (!is.finite(sd_pr) || sd_pr <= .Machine$double.eps) sd_pr <- 1

    .xform <- function(t, p) {
      cbind((t - mu_ta) / sd_ta, (p - mu_pr) / sd_pr)
    }
  } else {
    .xform <- function(t, p) cbind(t, p)
  }

  surf_xy <- .xform(surface_weights[[ta_col]], surface_weights[[pr_col]])
  if (any(!is.finite(surf_xy))) stop("surface_weights contains non-finite ta/pr values.", call. = FALSE)

  # ---- distance + mapping helpers ----
  .row_min_dist_and_argmin <- function(xy, surf_xy) {
    # returns list(min_d2, idx_min) for each row in xy
    n <- nrow(xy)
    min_d2 <- rep(Inf, n)
    idx_min <- rep(NA_integer_, n)

    for (s in seq.int(1L, n, by = chunk_size)) {
      e <- as.integer(min(s + chunk_size - 1L, n))
      block <- xy[s:e, , drop = FALSE]

      # brute force: compute squared distances to all surface points
      # memory: (block_size x M) can be large; do it in two passes to avoid full matrix
      # Here: loop surface points (O(block*M)) but minimal memory.
      bd2 <- rep(Inf, nrow(block))
      bidx <- rep(NA_integer_, nrow(block))

      for (j in seq_len(nrow(surf_xy))) {
        dx <- block[, 1] - surf_xy[j, 1]
        dy <- block[, 2] - surf_xy[j, 2]
        d2 <- dx * dx + dy * dy
        take <- d2 < bd2
        bd2[take] <- d2[take]
        bidx[take] <- j
      }

      min_d2[s:e] <- bd2
      idx_min[s:e] <- bidx
    }

    list(min_d = sqrt(min_d2), idx = idx_min)
  }

  .row_knn_weights <- function(xy, surf_xy, k, h) {
    # Returns list(idx = (n x k), w = (n x k), dist = (n x k))
    n <- nrow(xy)
    M <- nrow(surf_xy)
    if (k > M) stop("k cannot exceed number of surface points.", call. = FALSE)

    idx_mat <- matrix(NA_integer_, n, k)
    dist_mat <- matrix(NA_real_, n, k)

    for (i in seq_len(n)) {
      dx <- surf_xy[, 1] - xy[i, 1]
      dy <- surf_xy[, 2] - xy[i, 2]
      d2 <- dx * dx + dy * dy
      ord <- order(d2)[seq_len(k)]
      idx_mat[i, ] <- ord
      dist_mat[i, ] <- sqrt(d2[ord])
    }

    if (is.null(h)) {
      # median distance to kth neighbor across points
      h_use <- stats::median(dist_mat[, k])
      if (!is.finite(h_use) || h_use <= .Machine$double.eps) h_use <- 1
    } else {
      h_use <- h
      if (!is.finite(h_use) || h_use <= .Machine$double.eps) stop("kernel_bandwidth must be > 0.", call. = FALSE)
    }

    wker <- exp(-0.5 * (dist_mat / h_use)^2)
    # normalize kernel weights per row
    rs <- rowSums(wker)
    rs[rs <= .Machine$double.eps] <- NA_real_
    wker <- wker / rs

    list(idx = idx_mat, w = wker, dist = dist_mat, h = h_use)
  }

  # ---- score per group ----
  groups <- sort(unique(as.character(ensemble_ok[[group_col]])))
  # Only score groups that exist both in ensemble and surface columns
  groups <- intersect(groups, weight_cols)

  if (length(groups) == 0L) {
    stop("No overlapping groups between ensemble_data and surface_weights columns.", call. = FALSE)
  }

  res <- vector("list", length(groups))
  names(res) <- groups

  for (g in groups) {

    df_g <- ensemble_ok[as.character(ensemble_ok[[group_col]]) == g, , drop = FALSE]
    if (nrow(df_g) == 0L) next

    xy <- .xform(df_g[[ta_col]], df_g[[pr_col]])

    # scoring weights for aggregation
    if (is.null(weights_col)) {
      w_score <- rep(1, nrow(df_g))
    } else {
      w_score <- df_g[[weights_col]]
    }
    w_score <- w_score / sum(w_score)

    w_surface <- surface_weights[[g]]
    if (!is.numeric(w_surface) || length(w_surface) != nrow(surface_weights)) {
      stop("Surface weight column '", g, "' is not numeric or has wrong length.", call. = FALSE)
    }
    if (any(!is.finite(w_surface)) || any(w_surface < 0)) {
      stop("Surface weight column '", g, "' contains non-finite or negative values.", call. = FALSE)
    }

    # If not normalized, we still score but warn via flag
    sum_w <- sum(w_surface)
    normalized <- is.finite(sum_w) && abs(sum_w - 1) < 1e-6

    if (mapping == "nearest") {
      nn <- .row_min_dist_and_argmin(xy, surf_xy)
      p <- w_surface[nn$idx]
      p <- pmax(p, eps)

      mean_log_score <- sum(w_score * log(p))
      mean_nn_dist <- sum(w_score * nn$min_d)
      zero_mass_rate <- sum(w_score * (w_surface[nn$idx] <= eps))

      res[[g]] <- data.frame(
        group = g,
        mapping = "nearest",
        n = nrow(df_g),
        mean_log_score = mean_log_score,
        neg_log_score = -mean_log_score,
        mean_nn_distance = mean_nn_dist,
        zero_mass_rate = zero_mass_rate,
        surface_normalized = normalized,
        stringsAsFactors = FALSE
      )

    } else {

      knn <- .row_knn_weights(xy, surf_xy, k = k, h = kernel_bandwidth)
      # p_i = sum_j wker_ij * w_surface[idx_ij]
      p_raw <- rowSums(knn$w * matrix(w_surface[knn$idx], nrow = nrow(knn$idx)))
      p <- pmax(p_raw, eps)

      mean_log_score <- sum(w_score * log(p))
      mean_nn_dist <- sum(w_score * knn$dist[, 1])
      zero_mass_rate <- sum(w_score * (p_raw <= eps))

      res[[g]] <- data.frame(
        group = g,
        mapping = paste0("knn(k=", k, ",h=", format(knn$h, digits = 3), ")"),
        n = nrow(df_g),
        mean_log_score = mean_log_score,
        neg_log_score = -mean_log_score,
        mean_nn_distance = mean_nn_dist,
        zero_mass_rate = zero_mass_rate,
        surface_normalized = normalized,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, res)

  # ---- overall (weighted by group n) ----
  ok <- is.finite(out$mean_log_score)
  if (any(ok)) {
    w <- out$n[ok] / sum(out$n[ok])
    overall <- data.frame(
      group = "OVERALL",
      mapping = unique(out$mapping)[1],
      n = sum(out$n[ok]),
      mean_log_score = sum(w * out$mean_log_score[ok]),
      neg_log_score = -sum(w * out$mean_log_score[ok]),
      mean_nn_distance = sum(w * out$mean_nn_distance[ok]),
      zero_mass_rate = sum(w * out$zero_mass_rate[ok]),
      surface_normalized = all(out$surface_normalized[ok]),
      stringsAsFactors = FALSE
    )
    out <- rbind(out, overall)
  }

  out
}
