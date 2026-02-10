# ==============================================================================
# Extended Evaluation of Scenario Surface Weights
# ==============================================================================
#
# Goes beyond statistical fit to include decision-relevant metrics for
# climate stress testing and adaptation workflows.
#
# Evaluation dimensions:
#   1. Statistical calibration (log-scores)
#   2. Weight concentration (effective N, Gini, dominance)
#   3. Tail/stress coverage (user-defined regions, quantile-based)
#   4. Spatial coherence (centroid alignment, dispersion)
#   5. Robustness (stability under leave-one-out)
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Concentration Metrics
# ------------------------------------------------------------------------------

#' Effective number of scenarios (entropy-based)
#'
#' Measures how many scenarios effectively contribute to the distribution.
#' Range: 1 (all mass on one cell) to n (uniform).
#'
#' @param w Numeric vector of weights (will be normalized)
#' @return Effective N
#' @export
effective_n_entropy <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  w <- w / sum(w)
  # Shannon entropy -> effective N

  H <- -sum(w * log(w))
  exp(H)
}

#' Effective number of scenarios (HHI-based / Simpson)
#'
#' Inverse Herfindahl-Hirschman Index. Less sensitive to small weights than entropy.
#'
#' @param w Numeric vector of weights
#' @return Effective N
#' @export
effective_n_hhi <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  w <- w / sum(w)
  1 / sum(w^2)
}

#' Gini coefficient of weight distribution
#'
#' 0 = perfect equality (uniform), 1 = perfect inequality (all mass on one cell).
#'
#' @param w Numeric vector of weights
#' @return Gini coefficient
#' @export
weight_gini <- function(w) {
  w <- w[is.finite(w)]
  if (length(w) < 2) return(NA_real_)
  w <- sort(w)

  n <- length(w)
  # Gini formula
  numerator <- 2 * sum(seq_len(n) * w)
  denominator <- n * sum(w)
  gini <- (numerator / denominator) - (n + 1) / n
  gini
}

#' Maximum weight (dominance measure)
#'
#' High values indicate single-scenario dominance.
#'
#' @param w Numeric vector of weights
#' @return Maximum weight after normalization
#' @export
weight_max <- function(w) {
  w <- w[is.finite(w)]
  if (length(w) == 0) return(NA_real_)
  max(w) / sum(w)
}

#' Number of cells with non-negligible weight
#'
#' @param w Numeric vector of weights
#' @param threshold Minimum weight fraction to count as "active"
#' @return Count of active cells
#' @export
n_active_cells <- function(w, threshold = 0.01) {
  w <- w[is.finite(w)]
  if (length(w) == 0) return(NA_integer_)
  w <- w / sum(w)
  sum(w >= threshold)
}

# ------------------------------------------------------------------------------
# Tail/Stress Coverage Metrics
# ------------------------------------------------------------------------------
#' Weight in tail regions (quantile-based)
#'
#' Computes fraction of weight in joint distribution tails.
#'
#' @param w Numeric vector of weights
#' @param ta Grid temperature values
#' @param pr Grid precipitation values
#' @param tail_prob Quantile threshold (e.g., 0.1 = 10th/90th percentiles)
#' @return Named list with lower/upper tail weights
#' @export
tail_coverage_quantile <- function(w, ta, pr, tail_prob = 0.1) {

  if (length(w) != length(ta) || length(w) != length(pr)) {
    stop("w, ta, pr must have same length")
  }

  w <- w / sum(w[is.finite(w)])

  ta_lo <- quantile(ta, tail_prob, na.rm = TRUE)
  ta_hi <- quantile(ta, 1 - tail_prob, na.rm = TRUE)
  pr_lo <- quantile(pr, tail_prob, na.rm = TRUE)
  pr_hi <- quantile(pr, 1 - tail_prob, na.rm = TRUE)

  # Four corners
  hot_dry <- (ta >= ta_hi) & (pr <= pr_lo)
  hot_wet <- (ta >= ta_hi) & (pr >= pr_hi)
  cold_dry <- (ta <= ta_lo) & (pr <= pr_lo)
  cold_wet <- (ta <= ta_lo) & (pr >= pr_hi)

  # Any tail (union of margins)
  any_tail <- (ta <= ta_lo) | (ta >= ta_hi) | (pr <= pr_lo) | (pr >= pr_hi)

  list(
    hot_dry = sum(w[hot_dry], na.rm = TRUE),
    hot_wet = sum(w[hot_wet], na.rm = TRUE),
    cold_dry = sum(w[cold_dry], na.rm = TRUE),
    cold_wet = sum(w[cold_wet], na.rm = TRUE),
    any_tail = sum(w[any_tail], na.rm = TRUE),
    tail_quantile = tail_prob
  )
}

#' Weight in user-defined stress regions
#'
#' @param w Numeric vector of weights
#' @param ta Grid temperature values
#' @param pr Grid precipitation values
#' @param stress_regions Named list of regions, each with ta_min, ta_max, pr_min, pr_max
#' @return Named vector of weights per region
#' @export
stress_region_coverage <- function(w, ta, pr, stress_regions) {

  if (is.null(stress_regions) || length(stress_regions) == 0) {
    return(NULL)
  }

  w <- w / sum(w[is.finite(w)])

  out <- sapply(names(stress_regions), function(region_name) {
    region <- stress_regions[[region_name]]

    ta_min <- if (!is.null(region$ta_min)) region$ta_min else -Inf
    ta_max <- if (!is.null(region$ta_max)) region$ta_max else Inf
    pr_min <- if (!is.null(region$pr_min)) region$pr_min else -Inf
    pr_max <- if (!is.null(region$pr_max)) region$pr_max else Inf

    mask <- (ta >= ta_min) & (ta <= ta_max) & (pr >= pr_min) & (pr <= pr_max)
    sum(w[mask], na.rm = TRUE)
  })

  out
}

# ------------------------------------------------------------------------------
# Spatial Coherence Metrics
# ------------------------------------------------------------------------------

#' Weighted centroid of distribution
#'
#' @param w Weights
#' @param ta Temperature values
#' @param pr Precipitation values
#' @return Named vector c(ta = ..., pr = ...)
#' @keywords internal
weighted_centroid <- function(w, ta, pr) {
  w <- w / sum(w[is.finite(w)])
  c(ta = sum(w * ta, na.rm = TRUE), pr = sum(w * pr, na.rm = TRUE))
}

#' Data centroid (unweighted mean of observations)
#' @keywords internal
data_centroid <- function(ta_obs, pr_obs) {
  c(ta = mean(ta_obs, na.rm = TRUE), pr = mean(pr_obs, na.rm = TRUE))
}

#' Centroid shift: distance between weight centroid and data centroid
#'
#' Large values suggest the surface mass is offset from where data actually is.
#' Computed in standardized units to make T and P comparable.
#'
#' @param w Weights
#' @param ta_grid Grid temperature values
#' @param pr_grid Grid precipitation values
#' @param ta_obs Observed temperature values
#' @param pr_obs Observed precipitation values
#' @return Euclidean distance in standardized units
#' @export
centroid_shift <- function(w, ta_grid, pr_grid, ta_obs, pr_obs) {

  # Standardize using observation sd
  ta_sd <- sd(ta_obs, na.rm = TRUE)
  pr_sd <- sd(pr_obs, na.rm = TRUE)

  if (!is.finite(ta_sd) || ta_sd < .Machine$double.eps) ta_sd <- 1

  if (!is.finite(pr_sd) || pr_sd < .Machine$double.eps) pr_sd <- 1

  wc <- weighted_centroid(w, ta_grid / ta_sd, pr_grid / pr_sd)
  dc <- data_centroid(ta_obs / ta_sd, pr_obs / pr_sd)

  sqrt(sum((wc - dc)^2))
}

#' Weighted dispersion (spread of mass)
#'
#' Measures how spread out the weight distribution is around its centroid.
#' Higher = more spread, lower = concentrated.
#'
#' @param w Weights
#' @param ta Grid temperature values
#' @param pr Grid precipitation values
#' @return Weighted standard distance (standardized)
#' @export
weighted_dispersion <- function(w, ta, pr) {

  # Standardize
  ta_sd <- sd(ta, na.rm = TRUE)
  pr_sd <- sd(pr, na.rm = TRUE)
  if (!is.finite(ta_sd) || ta_sd < .Machine$double.eps) ta_sd <- 1
  if (!is.finite(pr_sd) || pr_sd < .Machine$double.eps) pr_sd <- 1

  ta_s <- ta / ta_sd
  pr_s <- pr / pr_sd

  w <- w / sum(w[is.finite(w)])

  wc <- weighted_centroid(w, ta_s, pr_s)

  # Weighted variance
  var_ta <- sum(w * (ta_s - wc["ta"])^2, na.rm = TRUE)
  var_pr <- sum(w * (pr_s - wc["pr"])^2, na.rm = TRUE)

  sqrt(var_ta + var_pr)
}

# ------------------------------------------------------------------------------
# Robustness Metrics
# ------------------------------------------------------------------------------

#' Leave-one-out weight stability
#'
#' Measures how much the weight surface changes when removing each observation.
#' High instability suggests overfitting to individual GCMs.
#'
#' Note: Requires a refit function. Returns NA if not provided.
#'
#' @param weight_fn Function that takes ensemble_data and returns weight surface
#' @param ensemble_data Full ensemble data
#' @param scenario_grid Evaluation grid
#' @param group Current group being evaluated
#' @param group_col Group column name
#' @return Mean absolute weight change across LOO iterations
#' @export
loo_stability <- function(weight_fn, ensemble_data, scenario_grid, group, group_col = "scenario") {

  if (is.null(weight_fn)) return(NA_real_)

  # Subset to group
  group_data <- ensemble_data[ensemble_data[[group_col]] == group, , drop = FALSE]
  n <- nrow(group_data)

  if (n < 5) return(NA_real_)  # Need enough data for LOO

  # Helper: extract weight vector from either a numeric vector (legacy)
  # or a long-format weight surface with a `weight` column.
  .extract_w <- function(x) {
    if (is.numeric(x) && !is.matrix(x) && !is.data.frame(x)) {
      return(as.numeric(x))
    }
    if (is.data.frame(x)) {
      if (!("weight" %in% names(x))) {
        stop("weight_fn returned a data.frame without a 'weight' column.", call. = FALSE)
      }
      if ("scenario" %in% names(x)) {
        x <- x[x$scenario == group, , drop = FALSE]
      }
      return(as.numeric(x$weight))
    }
    stop("weight_fn must return either a numeric weight vector or a data.frame with a 'weight' column.",
         call. = FALSE)
  }

  # Fit full model
  tryCatch({
    w_full <- .extract_w(weight_fn(group_data, scenario_grid))
    w_full <- w_full / sum(w_full, na.rm = TRUE)

    # LOO fits
    w_loo <- matrix(NA_real_, nrow = length(w_full), ncol = n)

    for (i in seq_len(n)) {
      w_i <- .extract_w(weight_fn(group_data[-i, , drop = FALSE], scenario_grid))
      w_loo[, i] <- w_i / sum(w_i, na.rm = TRUE)
    }

    # Mean absolute deviation from full fit
    mean(abs(w_loo - w_full), na.rm = TRUE)

  }, error = function(e) NA_real_)
}

# ------------------------------------------------------------------------------
# Main Extended Evaluation Function
# ------------------------------------------------------------------------------

#' Extended Evaluation of Scenario Surface Weights
#'
#' Computes statistical fit, concentration, tail coverage, spatial coherence,
#' and robustness metrics for climate scenario weight surfaces.
#'
#' @param weight_surface Data frame with grid coordinates and weight columns per group
#' @param ensemble_data Data frame with GCM projections
#' @param ta_col Temperature column name
#' @param pr_col Precipitation column name
#' @param group_col Group/scenario column name
#' @param mapping Mapping method for log-score: "nearest", "knn", or "bilinear"
#' @param k Number of neighbors for knn mapping
#' @param tail_quantile Quantile threshold for tail regions (default 0.1)
#' @param stress_regions Optional named list of user-defined stress regions
#' @param active_threshold Minimum weight to count as "active" cell
#' @param refit_fn Optional function for LOO stability (takes data, grid -> weights)
#' @param verbose Print progress?
#'
#' @return List with:
#'   - `by_group`: Data frame of per-group metrics
#'   - `overall`: Aggregated metrics across groups
#'   - `diagnostics`: Detailed diagnostic info
#'   - `recommendations`: Decision-support flags
#'
#' @export
evaluate_scenario_surface_weights <- function(
    weight_surface,
    ensemble_data,
    ta_col = "tavg",
    pr_col = "prcp",
    group_col = "scenario",
    mapping = c("knn", "nearest", "bilinear"),
    k = 5L,
    tail_quantile = 0.1,
    stress_regions = NULL,
    active_threshold = 0.01,
    refit_fn = NULL,
    verbose = TRUE
) {

  mapping <- match.arg(mapping)

  # ---------------------------------------------------------------------------
  # Validate long-format surface and align groups
  # ---------------------------------------------------------------------------

  required_surface_cols <- c(ta_col, pr_col, "scenario", "weight")
  missing_surface_cols <- setdiff(required_surface_cols, names(weight_surface))
  if (length(missing_surface_cols) > 0L) {
    stop(
      "weight_surface must be long-format with columns: ",
      paste(required_surface_cols, collapse = ", "),
      ". Missing: ", paste(missing_surface_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(weight_surface$weight)) {
    stop("weight_surface$weight must be numeric.", call. = FALSE)
  }
  if (any(!is.na(weight_surface$weight) & weight_surface$weight < 0)) {
    stop("weight_surface$weight must be non-negative.", call. = FALSE)
  }

  groups_in_data <- unique(ensemble_data[[group_col]])
  groups_in_surface <- unique(weight_surface$scenario)

  groups <- intersect(as.character(groups_in_data), as.character(groups_in_surface))
  if (length(groups) == 0L) {
    stop("No overlapping scenarios between ensemble_data[[group_col]] and weight_surface$scenario.",
         call. = FALSE)
  }

  # Establish a reference grid ordering using the first available group.
  .order_surface <- function(df) {
    df[order(df[[ta_col]], df[[pr_col]]), , drop = FALSE]
  }

  ref_grp <- groups[1]
  ref_df <- weight_surface[as.character(weight_surface$scenario) == ref_grp, , drop = FALSE]
  ref_df <- .order_surface(ref_df)

  if (nrow(ref_df) == 0L) {
    stop("weight_surface contains no rows for scenario '", ref_grp, "'.", call. = FALSE)
  }

  ta_grid <- ref_df[[ta_col]]
  pr_grid <- ref_df[[pr_col]]
  n_grid <- nrow(ref_df)

  # Scenario grid passed to refit_fn/LOO robustness: ta/pr grid only.
  scenario_grid_refit <- ref_df[, c(ta_col, pr_col), drop = FALSE]

  # ---------------------------------------------------------------------------
  # Per-group evaluation
  # ---------------------------------------------------------------------------

  results_list <- list()

  for (grp in groups) {

    if (isTRUE(verbose)) message("Evaluating group: ", grp)

    # Get weights and validate that the grid matches the reference grid
    df_w <- weight_surface[as.character(weight_surface$scenario) == grp, , drop = FALSE]
    df_w <- .order_surface(df_w)

    if (nrow(df_w) != n_grid) {
      stop("Scenario '", grp, "' has ", nrow(df_w), " grid rows; expected ", n_grid,
           " (must match across scenarios).", call. = FALSE)
    }

    if (!isTRUE(all.equal(df_w[[ta_col]], ta_grid, check.attributes = FALSE)) ||
        !isTRUE(all.equal(df_w[[pr_col]], pr_grid, check.attributes = FALSE))) {
      stop("Grid (", ta_col, ", ", pr_col, ") does not match across scenarios after ordering.",
           call. = FALSE)
    }

    w <- df_w$weight
    obs_data <- ensemble_data[as.character(ensemble_data[[group_col]]) == grp, , drop = FALSE]

    ta_obs <- obs_data[[ta_col]]
    pr_obs <- obs_data[[pr_col]]
    n_obs <- sum(is.finite(ta_obs) & is.finite(pr_obs))

    # --- Statistical Fit ---
    log_score <- tryCatch({
      compute_log_score(w, ta_grid, pr_grid, ta_obs, pr_obs, mapping = mapping, k = k)
    }, error = function(e) NA_real_)

    zero_mass <- tryCatch({
      compute_zero_mass_rate(w, ta_grid, pr_grid, ta_obs, pr_obs, mapping = mapping, k = k)
    }, error = function(e) NA_real_)

    # --- Concentration ---
    eff_n_ent <- effective_n_entropy(w)
    eff_n_hhi <- effective_n_hhi(w)
    gini <- weight_gini(w)
    max_w <- weight_max(w)
    n_active <- n_active_cells(w, threshold = active_threshold)

    # --- Tail Coverage ---
    tail_cov <- tail_coverage_quantile(w, ta_grid, pr_grid, tail_prob = tail_quantile)
    stress_cov <- stress_region_coverage(w, ta_grid, pr_grid, stress_regions)

    # --- Spatial Coherence ---
    c_shift <- centroid_shift(w, ta_grid, pr_grid, ta_obs, pr_obs)
    dispersion <- weighted_dispersion(w, ta_grid, pr_grid)

    # --- Robustness (if refit function provided) ---
    loo_stab <- if (!is.null(refit_fn)) {
      loo_stability(refit_fn, ensemble_data, scenario_grid_refit, grp, group_col)
    } else {
      NA_real_
    }

    # Compile results
    res <- data.frame(
      group = grp,
      n_obs = n_obs,
      n_grid = n_grid,

      # Statistical fit
      mean_log_score = log_score,
      zero_mass_rate = zero_mass,

      # Concentration
      effective_n_entropy = eff_n_ent,
      effective_n_hhi = eff_n_hhi,
      gini = gini,
      max_weight = max_w,
      n_active_cells = n_active,

      # Spatial
      centroid_shift = c_shift,
      dispersion = dispersion,

      # Robustness
      loo_stability = loo_stab,

      # Tail coverage (unpack)
      tail_hot_dry = tail_cov$hot_dry,
      tail_hot_wet = tail_cov$hot_wet,
      tail_cold_dry = tail_cov$cold_dry,
      tail_cold_wet = tail_cov$cold_wet,
      tail_any = tail_cov$any_tail,

      stringsAsFactors = FALSE
    )

    # Add stress regions if defined
    if (!is.null(stress_cov)) {
      for (region_name in names(stress_cov)) {
        res[[paste0("stress_", region_name)]] <- stress_cov[region_name]
      }
    }

    results_list[[grp]] <- res
  }

  by_group <- do.call(rbind, results_list)
  rownames(by_group) <- NULL

  # ---------------------------------------------------------------------------
  # Overall aggregation
  # ---------------------------------------------------------------------------

  overall <- data.frame(
    group = "OVERALL",
    n_obs = sum(by_group$n_obs),
    n_grid = n_grid,

    mean_log_score = mean(by_group$mean_log_score, na.rm = TRUE),
    zero_mass_rate = mean(by_group$zero_mass_rate, na.rm = TRUE),

    effective_n_entropy = mean(by_group$effective_n_entropy, na.rm = TRUE),
    effective_n_hhi = mean(by_group$effective_n_hhi, na.rm = TRUE),
    gini = mean(by_group$gini, na.rm = TRUE),
    max_weight = mean(by_group$max_weight, na.rm = TRUE),
    n_active_cells = mean(by_group$n_active_cells, na.rm = TRUE),

    centroid_shift = mean(by_group$centroid_shift, na.rm = TRUE),
    dispersion = mean(by_group$dispersion, na.rm = TRUE),

    loo_stability = mean(by_group$loo_stability, na.rm = TRUE),

    tail_hot_dry = mean(by_group$tail_hot_dry, na.rm = TRUE),
    tail_hot_wet = mean(by_group$tail_hot_wet, na.rm = TRUE),
    tail_cold_dry = mean(by_group$tail_cold_dry, na.rm = TRUE),
    tail_cold_wet = mean(by_group$tail_cold_wet, na.rm = TRUE),
    tail_any = mean(by_group$tail_any, na.rm = TRUE),

    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Decision-support recommendations
  # ---------------------------------------------------------------------------

  recommendations <- generate_recommendations(by_group, overall, n_grid)

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  structure(
    list(
      by_group = by_group,
      overall = overall,
      recommendations = recommendations,
      params = list(
        mapping = mapping,
        k = k,
        tail_quantile = tail_quantile,
        stress_regions = stress_regions,
        active_threshold = active_threshold
      )
    ),
    class = c("scenario_weight_evaluation", "list")
  )
}

# ------------------------------------------------------------------------------
# Recommendation Generator
# ------------------------------------------------------------------------------

#' Generate decision-support recommendations from evaluation metrics
#' @keywords internal
generate_recommendations <- function(by_group, overall, n_grid) {

  flags <- list()

  # --- Concentration flags ---

  eff_n <- overall$effective_n_entropy
  if (!is.na(eff_n)) {
    if (eff_n < 3) {
      flags$concentration <- list(
        level = "warning",
        message = sprintf(
          "Very high concentration: effective N = %.1f (only ~%.0f scenarios matter). Consider regularization or flatter method.",
          eff_n, eff_n
        )
      )
    } else if (eff_n < 5) {
      flags$concentration <- list(
        level = "note",
        message = sprintf(
          "Moderate concentration: effective N = %.1f. Verify conclusions aren't driven by 1-2 scenarios.",
          eff_n
        )
      )
    }
  }

  # Max weight dominance
  max_w <- overall$max_weight
  if (!is.na(max_w) && max_w > 0.5) {
    flags$dominance <- list(
      level = "warning",
      message = sprintf(
        "Single-cell dominance: max weight = %.1f%%. Results may be fragile.",
        max_w * 100
      )
    )
  }

  # --- Tail coverage flags ---

  tail_any <- overall$tail_any
  if (!is.na(tail_any)) {
    # Expected tail coverage under uniform would be ~36% for 10% tails on each margin
    # (1 - 0.8^2 for any margin in tail)
    if (tail_any < 0.05) {
      flags$tail_coverage <- list(
        level = "warning",
        message = sprintf(
          "Negligible tail coverage (%.1f%%). Extreme scenarios effectively ignored. Consider tail-aware method.",
          tail_any * 100
        )
      )
    } else if (tail_any < 0.15) {
      flags$tail_coverage <- list(
        level = "note",
        message = sprintf(
          "Low tail coverage (%.1f%%). Verify extreme scenarios are appropriately weighted for risk assessment.",
          tail_any * 100
        )
      )
    }
  }

  # Asymmetric tail coverage (potential bias)
  hot_dry <- overall$tail_hot_dry
  cold_wet <- overall$tail_cold_wet
  if (!is.na(hot_dry) && !is.na(cold_wet)) {
    ratio <- if (cold_wet > 0.001) hot_dry / cold_wet else Inf
    if (is.finite(ratio) && (ratio > 10 || ratio < 0.1)) {
      flags$tail_asymmetry <- list(
        level = "note",
        message = sprintf(
          "Asymmetric tail weights: hot-dry = %.2f%%, cold-wet = %.2f%%. Verify this reflects physics, not method artifact.",
          hot_dry * 100, cold_wet * 100
        )
      )
    }
  }

  # --- Spatial coherence flags ---

  c_shift <- overall$centroid_shift
  if (!is.na(c_shift) && c_shift > 0.5) {
    flags$centroid <- list(
      level = "note",
      message = sprintf(
        "Centroid shift = %.2f SD. Weight distribution offset from data center.",
        c_shift
      )
    )
  }

  # --- Sample size warning ---

  min_n <- min(by_group$n_obs, na.rm = TRUE)
  if (min_n < 20) {
    flags$sample_size <- list(
      level = "note",
      message = sprintf(
        "Small sample sizes (min n = %d). Parametric methods (MVN) may be more reliable than KDE/copula.",
        min_n
      )
    )
  }

  # --- Active cells ---

  n_active <- overall$n_active_cells
  if (!is.na(n_active) && n_active < 3) {
    flags$sparsity <- list(
      level = "warning",
      message = sprintf(
        "Only %.0f cells have weight > 1%%. Surface is nearly degenerate.",
        n_active
      )
    )
  }

  # --- Summary verdict ---

  n_warnings <- sum(sapply(flags, function(f) f$level == "warning"))
  n_notes <- sum(sapply(flags, function(f) f$level == "note"))

  if (n_warnings > 0) {
    verdict <- "REVIEW: Method may not be suitable for robust stress testing"
  } else if (n_notes > 1) {
    verdict <- "ACCEPTABLE with caveats: Review flagged issues"
  } else {
    verdict <- "GOOD: No major concerns identified"
  }

  list(
    flags = flags,
    n_warnings = n_warnings,
    n_notes = n_notes,
    verdict = verdict
  )
}

# ------------------------------------------------------------------------------
# Helper: Log-score computation (simplified, assumes external implementation)
# ------------------------------------------------------------------------------

#' Compute log-score via KNN mapping
#' @keywords internal
compute_log_score <- function(w, ta_grid, pr_grid, ta_obs, pr_obs, mapping = "knn", k = 5) {

  w <- w / sum(w, na.rm = TRUE)

  n_obs <- length(ta_obs)
  log_scores <- numeric(n_obs)

  for (i in seq_len(n_obs)) {
    # Distance to all grid points
    d <- sqrt((ta_grid - ta_obs[i])^2 + (pr_grid - pr_obs[i])^2)

    if (mapping == "nearest") {
      idx <- which.min(d)
      p <- w[idx]
    } else if (mapping == "knn") {
      idx <- order(d)[1:min(k, length(d))]
      # Inverse distance weights
      d_k <- d[idx]
      d_k[d_k < 1e-10] <- 1e-10
      idw <- 1 / d_k
      idw <- idw / sum(idw)
      p <- sum(w[idx] * idw)
    } else {
      stop("Unsupported mapping: ", mapping)
    }

    log_scores[i] <- if (p > 0) log(p) else -Inf
  }

  # Handle -Inf gracefully
  log_scores[!is.finite(log_scores)] <- -20  # floor
  mean(log_scores)
}

#' Compute zero-mass rate
#' @keywords internal
compute_zero_mass_rate <- function(w, ta_grid, pr_grid, ta_obs, pr_obs, mapping = "knn", k = 5) {

  w <- w / sum(w, na.rm = TRUE)
  n_obs <- length(ta_obs)
  zero_count <- 0

  for (i in seq_len(n_obs)) {
    d <- sqrt((ta_grid - ta_obs[i])^2 + (pr_grid - pr_obs[i])^2)

    if (mapping == "nearest") {
      idx <- which.min(d)
      p <- w[idx]
    } else if (mapping == "knn") {
      idx <- order(d)[1:min(k, length(d))]
      d_k <- d[idx]
      d_k[d_k < 1e-10] <- 1e-10
      idw <- 1 / d_k
      idw <- idw / sum(idw)
      p <- sum(w[idx] * idw)
    }

    if (p < 1e-10) zero_count <- zero_count + 1
  }

  zero_count / n_obs
}

# ------------------------------------------------------------------------------
# Print method
# ------------------------------------------------------------------------------

#' @export
print.scenario_weight_evaluation <- function(x, ...) {

  cat("=== Scenario Weight Evaluation ===\n\n")

  cat("-- Summary --\n")
  cat(sprintf("Groups: %d | Grid cells: %d | Total obs: %d\n",
              nrow(x$by_group), x$by_group$n_grid[1], x$overall$n_obs))
  cat("\n")

  cat("-- Key Metrics (Overall) --\n")
  cat(sprintf("  Log-score:       %.3f  (higher = better fit)\n", x$overall$mean_log_score))
  cat(sprintf("  Effective N:     %.1f   (entropy) / %.1f (HHI)\n",
              x$overall$effective_n_entropy, x$overall$effective_n_hhi))
  cat(sprintf("  Gini:            %.3f  (0=uniform, 1=concentrated)\n", x$overall$gini))
  cat(sprintf("  Max weight:      %.1f%%\n", x$overall$max_weight * 100))
  cat(sprintf("  Active cells:    %.0f\n", x$overall$n_active_cells))
  cat(sprintf("  Tail coverage:   %.1f%%\n", x$overall$tail_any * 100))
  cat(sprintf("  Centroid shift:  %.2f SD\n", x$overall$centroid_shift))
  cat("\n")

  cat("-- Recommendations --\n")
  cat(sprintf("Verdict: %s\n", x$recommendations$verdict))
  if (length(x$recommendations$flags) > 0) {
    for (name in names(x$recommendations$flags)) {
      flag <- x$recommendations$flags[[name]]
      prefix <- if (flag$level == "warning") "[!]" else "[i]"
      cat(sprintf("  %s %s\n", prefix, flag$message))
    }
  }
  cat("\n")

  invisible(x)
}

# ------------------------------------------------------------------------------
# Comparison helper
# ------------------------------------------------------------------------------

#' Compare multiple weighting methods
#'
#' @param ... Named evaluation results from evaluate_scenario_surface_weights()
#' @return Comparison data frame and recommendation
#' @export
compare_methods <- function(...) {

  evals <- list(...)
  if (is.null(names(evals)) || any(names(evals) == "")) {
    stop("All arguments must be named (e.g., kde = eval_kde, mvn = eval_mvn)")
  }

  # Extract overall metrics
  comparison <- do.call(rbind, lapply(names(evals), function(method) {
    ov <- evals[[method]]$overall
    data.frame(
      method = method,
      log_score = ov$mean_log_score,
      eff_n_entropy = ov$effective_n_entropy,
      eff_n_hhi = ov$effective_n_hhi,
      gini = ov$gini,
      max_weight = ov$max_weight,
      n_active = ov$n_active_cells,
      tail_coverage = ov$tail_any,
      centroid_shift = ov$centroid_shift,
      dispersion = ov$dispersion,
      n_warnings = evals[[method]]$recommendations$n_warnings,
      stringsAsFactors = FALSE
    )
  }))

  # Rank methods on key criteria
  comparison$rank_logscore <- rank(-comparison$log_score)  # higher is better
  comparison$rank_effn <- rank(-comparison$eff_n_entropy)  # higher is better
  comparison$rank_tail <- rank(-comparison$tail_coverage)  # higher is better
  comparison$rank_warnings <- rank(comparison$n_warnings)  # lower is better

  # Simple aggregate rank
  comparison$rank_overall <- rowMeans(comparison[, c("rank_logscore", "rank_effn", "rank_tail", "rank_warnings")])

  # Sort by overall rank
  comparison <- comparison[order(comparison$rank_overall), ]

  structure(
    list(
      comparison = comparison,
      best_statistical = comparison$method[which.min(comparison$rank_logscore)],
      best_robust = comparison$method[which.min(comparison$rank_effn)],
      best_overall = comparison$method[1]
    ),
    class = c("method_comparison", "list")
  )
}

#' @export
print.method_comparison <- function(x, ...) {

  cat("=== Method Comparison ===\n\n")

  # Print key columns
  cols <- c("method", "log_score", "eff_n_entropy", "gini", "tail_coverage", "n_warnings", "rank_overall")
  print(x$comparison[, cols], row.names = FALSE)

  cat("\n")
  cat(sprintf("Best statistical fit:  %s\n", x$best_statistical))
  cat(sprintf("Best robustness:       %s\n", x$best_robust))
  cat(sprintf("Best overall:          %s\n", x$best_overall))
  cat("\n")

  invisible(x)
}
