#' Calculate Climate-Informed Scenario Weights (KDE; auto bandwidth; optional family weighting)
#'
#' @description
#' For each group (e.g., SSP), estimates a smooth 2D kernel density surface over
#' (temperature change, precipitation change) using Gaussian KDE, evaluates it
#' on a user-provided stress-test grid, and returns normalized weights.
#'
#' If `bw = NULL`, bandwidths are chosen automatically per group using a
#' grid-aware rule:
#' \deqn{bw = max(k * grid_step, alpha * bw_nrd)}
#' where `grid_step` is the median spacing in `scenario_grid` for each dimension,
#' `bw_nrd` is a data-driven bandwidth estimate for each dimension, `k` controls
#' smoothing relative to grid spacing, and `alpha` scales the data-driven bandwidth.
#'
#' Optionally, observations can be downweighted by model family (near-duplicate control)
#' so that each family contributes equal total influence within each group.
#'
#' @param ensemble_data Data frame with temperature, precipitation, and group columns.
#' @param scenario_grid Data frame with all stress-test combinations (must include `ta_col`, `pr_col`).
#' @param pr_col Character. Precipitation column name. Default `"prcp"`.
#' @param ta_col Character. Temperature column name. Default `"tavg"`.
#' @param group_col Character. Grouping column name (e.g., `"scenario"`). Default `"scenario"`.
#' @param bw Numeric length-2 bandwidth vector `c(bw_ta, bw_pr)`, or `NULL` (default) for auto.
#' @param k Numeric length-2. Grid-step multipliers for auto bw. Default `c(1.5, 2.0)`.
#' @param alpha Numeric scalar. Multiplier on data-driven bandwidth in auto bw. Default `1.0`.
#' @param bw_min Numeric length-2. Minimum allowed bandwidths. Default `c(0, 0)`.
#' @param bw_max Numeric length-2. Maximum allowed bandwidths. Default `c(Inf, Inf)`.
#' @param min_samples Integer. Minimum complete observations per group. Default `5`.
#' @param normalize Logical. If TRUE (default), normalize weights within each group to sum to 1.
#' @param use_family_weights Logical. If TRUE, downweight near-duplicates by `family_col`.
#' @param family_col Character. Column name giving model family. Default `"model_family"`.
#' @param chunk_size Integer. Grid chunk size for memory control. Default `5000`.
#'
#' @return
#' Data frame identical to `scenario_grid` plus one weight column per group (alphabetically ordered).
#' Attribute `"skipped_groups"` lists groups skipped due to insufficient data or degeneracy.
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
#' w <- estimate_scenario_probabilities_kde(
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
#' # --- Example 2: Fixed user bandwidth (overrides auto for all groups) ---
#' w_fixed <- estimate_scenario_probabilities_kde(
#'   ensemble_data = ensemble_data,
#'   scenario_grid = scenario_grid,
#'   bw            = c(0.35, 0.90),
#'   verbose       = FALSE
#' )
#' colSums(w_fixed[, c("SSP1", "SSP2")])
#'
#' # --- Example 3: Family downweighting to reduce near-duplicate influence ---
#' ensemble_data$model_family <- rep(c("A", "A", "B"), length.out = nrow(ensemble_data))
#'
#' w_family <- estimate_scenario_probabilities_kde(
#'   ensemble_data        = ensemble_data,
#'   scenario_grid        = scenario_grid,
#'   use_family_weights   = TRUE,
#'   family_col           = "model_family",
#'   verbose              = FALSE
#' )
#' colSums(w_family[, c("SSP1", "SSP2")])
#'
#' # --- Example 4: Groups can be skipped (insufficient data / degenerate variance) ---
#' ensemble_small <- rbind(
#'   ensemble_data,
#'   data.frame(scenario = "TOO_SMALL", tavg = c(1, 2), prcp = c(0, 1)),
#'   data.frame(scenario = "DEGENERATE", tavg = rep(2, 10), prcp = rnorm(10))
#' )
#'
#' w_skip <- estimate_scenario_probabilities_kde(
#'   ensemble_data = ensemble_small,
#'   scenario_grid = scenario_grid,
#'   min_samples   = 5L,
#'   verbose       = FALSE
#' )
#' attr(w_skip, "skipped_groups")
#'
#' # --- Example 5: Group names are made syntactically valid in output columns ---
#' ensemble_named <- ensemble_data
#' ensemble_named$scenario <- ifelse(ensemble_named$scenario == "SSP1", "SSP2-4.5", "SSP5-8.5")
#' w_names <- estimate_scenario_probabilities_kde(
#'   ensemble_data = ensemble_named,
#'   scenario_grid = scenario_grid,
#'   verbose       = FALSE
#' )
#' names(w_names)[-(1:2)]
#'
#' @importFrom MASS bandwidth.nrd
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
estimate_scenario_probabilities_kde <- function(
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
    use_family_weights = FALSE,
    family_col = "model_family",
    chunk_size = 5000L,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  .grid_step <- function(x) {
    x <- sort(unique(x))
    if (length(x) < 2) return(NA_real_)
    dx <- diff(x)
    dx <- dx[is.finite(dx) & dx > .Machine$double.eps]
    if (length(dx) == 0) return(NA_real_)
    stats::median(dx)
  }

  .normalize_vec <- function(x) {
    s <- sum(x)
    if (is.finite(s) && s > .Machine$double.eps) x / s else x
  }

  .compute_family_weights <- function(family_id) {
    family_id <- as.character(family_id)
    family_id[is.na(family_id)] <- "NA_FAMILY"
    fam_tab <- table(family_id)
    n_fam <- length(fam_tab)
    w_fam <- 1 / n_fam
    w_i <- w_fam / as.numeric(fam_tab[family_id]) # split within family
    as.numeric(w_i)
  }

  .validate_bw_vec <- function(bw_vec) {
    is.numeric(bw_vec) && length(bw_vec) == 2 &&
      all(is.finite(bw_vec)) && all(bw_vec > 0)
  }

  .pick_bw_auto <- function(ta, pr, dT, dP, k, alpha, bw_min, bw_max) {
    bw_nrd <- tryCatch(
      c(MASS::bandwidth.nrd(ta), MASS::bandwidth.nrd(pr)),
      error = function(e) c(NA_real_, NA_real_)
    )

    # grid-aware floor
    bw_grid <- c(k[1] * dT, k[2] * dP)

    # if grid step is NA (degenerate grid), fall back to nrd only
    if (!is.finite(bw_grid[1])) bw_grid[1] <- 0
    if (!is.finite(bw_grid[2])) bw_grid[2] <- 0

    # combine
    bw_use <- pmax(bw_grid, alpha * bw_nrd, na.rm = TRUE)

    # clamp
    bw_use <- pmax(bw_use, bw_min)
    bw_use <- pmin(bw_use, bw_max)

    # final sanity
    if (!.validate_bw_vec(bw_use)) return(NULL)
    bw_use
  }

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

  if (!is.logical(normalize) || length(normalize) != 1) {
    stop("normalize must be a single logical.", call. = FALSE)
  }

  if (!is.numeric(min_samples) || length(min_samples) != 1 || min_samples < 3) {
    stop("min_samples must be a single integer >= 3.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  if (!is.null(bw) && !.validate_bw_vec(bw)) {
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
  if (any(bw_max <= 0)) {
    stop("bw_max must be > 0 (use Inf for no upper bound).", call. = FALSE)
  }


  if (any(bw_max <= bw_min)) {
    stop("bw_max must be greater than bw_min for both dimensions.", call. = FALSE)
  }

  if (use_family_weights && !(family_col %in% names(ensemble_data))) {
    stop("use_family_weights=TRUE but family_col not found in ensemble_data: ", family_col,
         call. = FALSE)
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 100) {
    stop("chunk_size must be a single integer >= 100.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  # ---------------------------------------------------------------------------
  # Precompute grid step (global; applies to all groups consistently)
  # ---------------------------------------------------------------------------

  dT <- .grid_step(scenario_grid[[ta_col]])
  dP <- .grid_step(scenario_grid[[pr_col]])

  # ---------------------------------------------------------------------------
  # Prepare output and loop
  # ---------------------------------------------------------------------------

  out <- scenario_grid
  groups <- sort(unique(as.character(ensemble_data[[group_col]])))
  skipped <- character(0)

  grid_ta <- scenario_grid[[ta_col]]
  grid_pr <- scenario_grid[[pr_col]]
  n_grid  <- nrow(scenario_grid)

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

    # Bandwidth choice
    bw_use <- bw
    if (is.null(bw_use)) {
      bw_use <- .pick_bw_auto(
        ta = ta, pr = pr,
        dT = dT, dP = dP,
        k = k, alpha = alpha,
        bw_min = bw_min, bw_max = bw_max
      )

      if (is.null(bw_use)) {
        warning("Skipping group '", grp, "' (auto bandwidth selection failed).", call. = FALSE)
        skipped <- c(skipped, grp)
        next
      }
    }

    if (verbose) {
      msg <- if (is.null(bw)) "auto" else "user"
      message(
        sprintf(
          "[KDE] Group=%s | bw_%s = (ta=%.3f, pr=%.3f)",
          grp, msg, bw_use[1], bw_use[2]
        )
      )
    }

    bw_ta <- bw_use[1]
    bw_pr <- bw_use[2]

    # Observation weights (uniform or family-weighted)
    if (use_family_weights) {
      fam <- df[[family_col]][ok]
      w_obs <- .compute_family_weights(fam)
      w_obs <- w_obs / sum(w_obs)  # sum-to-one across obs
    } else {
      w_obs <- rep(1 / length(ta), length(ta))    # sum-to-one across obs
    }

    # KDE evaluation (weighted Gaussian product kernel)
    n_obs <- length(ta)
    dens <- numeric(n_grid)

    inv_bw_ta <- 1 / bw_ta
    inv_bw_pr <- 1 / bw_pr
    norm_const <- inv_bw_ta * inv_bw_pr

    idx_starts <- seq.int(1L, n_grid, by = chunk_size)

    for (s in idx_starts) {
      e <- as.integer(min(s + chunk_size - 1L, n_grid))

      g_ta <- grid_ta[s:e]
      g_pr <- grid_pr[s:e]

      acc <- numeric(length(g_ta))
      for (i in seq_len(n_obs)) {
        u <- (g_ta - ta[i]) * inv_bw_ta
        v <- (g_pr - pr[i]) * inv_bw_pr
        acc <- acc + w_obs[i] * stats::dnorm(u) * stats::dnorm(v)
      }
      dens[s:e] <- acc * norm_const
    }

    if (normalize) dens <- .normalize_vec(dens)

    out[[make.names(grp)]] <- dens
  }

  # Ensure deterministic column order: grid cols first, then alphabetic group cols
  base_cols <- names(scenario_grid)
  weight_cols <- setdiff(names(out), base_cols)
  out <- out[, c(base_cols, sort(weight_cols)), drop = FALSE]

  if (length(skipped) > 0) {
    message("Skipped ", length(skipped), " group(s): ", paste(skipped, collapse = ", "))
    attr(out, "skipped_groups") <- skipped
  } else {
    attr(out, "skipped_groups") <- character(0)
  }

  out
}


################################################################################
################################################################################
################################################################################


#' Calculate Climate-Informed Scenario Weights
#'
#' @description
#' Fits a bivariate normal distribution to mean temperature and precipitation
#' changes from an ensemble of GCM projections, then evaluates normalized
#' probability densities for a grid of stress-test scenarios.
#' These densities serve as climate-informed weighting factors that can be used
#' to probabilistically weight scenario outcomes in stress-testing or
#' adaptation pathway analyses.
#'
#' @param ensemble_data A data frame containing GCM ensemble results with
#'   temperature, precipitation, and a grouping column (e.g., `scenario` or
#'   `model`).
#' @param scenario_grid A data frame (typically created with
#'   [tidyr::expand_grid()]) containing combinations of temperature and
#'   precipitation change values for which density weights are computed.
#' @param pr_col Character string. Name of the precipitation column in
#'   `ensemble_data` and `scenario_grid`. Default is `"prcp"`.
#' @param ta_col Character string. Name of the temperature column in
#'   `ensemble_data` and `scenario_grid`. Default is `"tavg"`.
#' @param group_col Character string. Column name in `ensemble_data` defining
#'   groups for which separate bivariate distributions will be fitted
#'   (typically `"scenario"` or `"model"`). Default is `"scenario"`.
#' @param robust Logical; if `TRUE` (default), use a robust covariance estimate
#'   via [MASS::cov.trob()]. If `FALSE`, use the classical covariance with
#'   [stats::cov()]. Robust estimation requires at least 5 observations per
#'   group; classical estimation requires at least 3.
#' @param normalize Logical; if `TRUE` (default), densities for each group are
#'   normalized so they sum to one, creating proper probability weights.
#'
#' @return
#' A data frame identical to `scenario_grid`, with one additional column per
#' group in `ensemble_data`. Each added column contains the normalized bivariate
#' density values representing the relative likelihood (weight) of each
#' stress-test combination under that group's fitted climate distribution.
#'
#' The returned object includes a `"skipped_groups"` attribute listing any
#' groups that were excluded due to insufficient data or estimation issues.
#'
#' @details
#' For each group in `ensemble_data`, the function:
#' \enumerate{
#'   \item Extracts temperature and precipitation change values.
#'   \item Estimates the mean vector and covariance matrix (robust or classical).
#'   \item Checks for singular or degenerate covariance matrices.
#'   \item Evaluates the bivariate normal density for each point in
#'         `scenario_grid` using [mvtnorm::dmvnorm()].
#'   \item Optionally normalizes densities so that their sum equals one.
#' }
#'
#' The output is particularly useful when combining climate projections with
#' scenario-based stress tests — for example, when weighting hydrological or
#' economic model outputs by the relative likelihood of each climate change
#' condition.
#'
#' **Robust vs. Classical Estimation:**
#' Robust covariance estimation (via `MASS::cov.trob()`) is less sensitive to
#' outliers and is recommended when the ensemble contains potential outliers.
#' Classical estimation uses standard `stats::cov()` and may be preferable for
#' small, well-behaved ensembles.
#'
#' @note
#' Groups will be skipped with a warning if they have:
#' \itemize{
#'   \item Insufficient observations (< 5 for robust, < 3 for classical)
#'   \item Singular or near-singular covariance matrices
#'   \item Zero variance in temperature or precipitation
#'   \item Covariance estimation or density computation errors
#' }
#' The returned data frame will not include weight columns for skipped groups.
#' A message will report the number and names of skipped groups, and they are
#' also available via `attr(result, "skipped_groups")`.
#'
#' Column names in the output are derived from group values and sanitized using
#' [base::make.names()] to ensure valid R identifiers.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(tidyr)
#'
#' # Load GCM ensemble data
#' gcm_data <- readr::read_csv("Rhine/data/CMIP6/annual_change_scalar_stats_summary_mean.csv") %>%
#'   filter(horizon == "near") %>%
#'   rename(prcp = precip, tavg = temp)
#'
#' # Create stress-test grid
#' stress_grid <- expand_grid(
#'   tavg = seq(0, 6, 1),
#'   prcp = seq(-30, 30, 10)
#' )
#'
#' # Calculate weights with robust estimation
#' weights <- estimate_scenario_probabilities(
#'   ensemble_data = gcm_data,
#'   scenario_grid = stress_grid,
#'   group_col = "scenario",
#'   robust = TRUE
#' )
#'
#' # Check for skipped groups
#' attr(weights, "skipped_groups")
#'
#' # Use weights for scenario analysis
#' weighted_outcomes <- stress_grid %>%
#'   left_join(model_outputs, by = c("tavg", "prcp")) %>%
#'   mutate(
#'     weighted_impact = impact * weights$SSP245,
#'     across(starts_with("SSP"), ~ . * impact, .names = "weighted_{.col}")
#'   )
#' }
#'
#' @importFrom MASS cov.trob
#' @importFrom mvtnorm dmvnorm
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
#'
estimate_scenario_probabilities <- function(
    ensemble_data,
    scenario_grid,
    pr_col = "prcp",
    ta_col = "tavg",
    group_col = "scenario",
    robust = TRUE,
    normalize = TRUE) {

  # --- Input validation ---
  required_cols <- c(ta_col, pr_col, group_col)
  missing_cols <- setdiff(required_cols, names(ensemble_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ensemble_data: ",
         paste(missing_cols, collapse = ", "))
  }

  if (!all(c(ta_col, pr_col) %in% names(scenario_grid))) {
    stop("Both temperature and precipitation columns must exist in scenario_grid.")
  }

  if (nrow(scenario_grid) == 0) {
    stop("scenario_grid cannot be empty.")
  }

  if (!is.logical(robust) || !is.logical(normalize)) {
    stop("Arguments 'robust' and 'normalize' must be logical.")
  }

  # --- Prepare output ---
  out <- scenario_grid
  groups <- unique(ensemble_data[[group_col]])
  skipped <- character(0)

  # Minimum samples required
  min_samples <- if (robust) 5 else 3

  # --- Loop over groups ---
  for (grp in groups) {

    # Filter data for current group
    df <- dplyr::filter(ensemble_data, .data[[group_col]] == grp)

    # Check for sufficient observations
    if (nrow(df) < min_samples) {
      warning("Skipping group '", grp, "' (need at least ", min_samples,
              " observations for ", if (robust) "robust" else "classical",
              " estimation).")
      skipped <- c(skipped, as.character(grp))
      next
    }

    # Extract relevant columns
    df_subset <- df[, c(ta_col, pr_col)]

    # --- Covariance estimation ---
    result <- tryCatch({
      if (robust) {
        fit <- MASS::cov.trob(df_subset)
        list(mu = fit$center, Sigma = fit$cov)
      } else {
        # Remove NA values for classical estimation
        df_clean <- na.omit(df_subset)
        if (nrow(df_clean) < min_samples) {
          stop("Insufficient non-NA observations after removing missing values")
        }
        list(
          mu = colMeans(df_clean),
          Sigma = stats::cov(df_clean)
        )
      }
    }, error = function(e) {
      warning("Covariance estimation failed for group '", grp, "': ", e$message)
      NULL
    })

    # Skip if estimation failed
    if (is.null(result)) {
      skipped <- c(skipped, as.character(grp))
      next
    }

    mu <- result$mu
    Sigma <- result$Sigma

    # --- Validate covariance matrix ---
    # Check for NA values
    if (any(is.na(Sigma))) {
      warning("Skipping group '", grp, "' (covariance matrix contains NA values).")
      skipped <- c(skipped, as.character(grp))
      next
    }

    # Check for zero or negative variance
    if (any(diag(Sigma) < .Machine$double.eps)) {
      warning("Skipping group '", grp,
              "' (zero or negative variance in one or more variables).")
      skipped <- c(skipped, as.character(grp))
      next
    }

    # Check for singularity
    det_Sigma <- tryCatch(det(Sigma), error = function(e) NA)
    if (is.na(det_Sigma) || det_Sigma <= .Machine$double.eps) {
      warning("Skipping group '", grp,
              "' (singular or near-singular covariance matrix).")
      skipped <- c(skipped, as.character(grp))
      next
    }

    # --- Compute bivariate normal densities ---
    dens <- tryCatch(
      mvtnorm::dmvnorm(
        scenario_grid[, c(ta_col, pr_col)],
        mean = mu,
        sigma = Sigma
      ),
      error = function(e) {
        warning("Density computation failed for group '", grp, "': ", e$message)
        NULL
      }
    )

    # Skip if density computation failed
    if (is.null(dens)) {
      skipped <- c(skipped, as.character(grp))
      next
    }

    # --- Normalize densities ---
    if (normalize) {
      dens_sum <- sum(dens)
      if (dens_sum > .Machine$double.eps) {
        dens <- dens / dens_sum
      } else {
        warning("Zero total density for group '", grp,
                "'. Using unnormalized values.")
      }
    }

    # --- Add to output with safe column name ---
    col_name <- make.names(as.character(grp))
    out[[col_name]] <- dens
  }

  # --- Report skipped groups ---
  if (length(skipped) > 0) {
    message("Skipped ", length(skipped), " group(s): ",
            paste(skipped, collapse = ", "))
    attr(out, "skipped_groups") <- skipped
  } else {
    attr(out, "skipped_groups") <- character(0)
  }

  return(out)
}




estimate_grid_probabilities <- function(
    gcm_df_w,
    stress_test_grid,
    scenario_col = "scenario",
    ta_col = "tavg",
    pr_col = "prcp",
    weight_col = "w_genealogy",
    ...) {

  prob_wide <- estimate_scenario_probabilities_kde(
    ensemble_data = gcm_df_w,
    scenario_grid = stress_test_grid[, c(ta_col, pr_col)],
    weight_col = weight_col,
    group_col = scenario_col,
    ...
  )

  # convert to long: one row per grid × scenario
  scenarios <- setdiff(names(prob_wide), c(ta_col, pr_col))

  out <- do.call(
    rbind,
    lapply(scenarios, function(sc) {
      data.frame(
        scenario = sc,
        stress_test_grid[, c(ta_col, pr_col)],
        prob = prob_wide[[sc]],
        stringsAsFactors = FALSE
      )
    })
  )

  out
}



# ==============================================================================
# Weighted ensemble reporting utilities (single script)
# - Weighted quantiles
# - Effective ensemble size
# - Exceedance probabilities
# - Standard reporting table per scenario (e.g., SSP)
# ==============================================================================


# ------------------------------------------------------------------------------
# Weighted quantile (empirical, step-function)
# ------------------------------------------------------------------------------
# Returns weighted quantiles for numeric vector x with non-negative weights w.
# Notes:
# - Works with any weights (need not already sum to 1).
# - Implements the "smallest x where CDF >= p" convention.
weighted_quantile <- function(x, w, probs = c(0.10, 0.25, 0.50, 0.75, 0.90)) {
  if (!is.numeric(x) || !is.numeric(w)) {
    stop("x and w must be numeric.", call. = FALSE)
  }
  if (!is.numeric(probs) || any(probs < 0 | probs > 1)) {
    stop("probs must be numeric in [0, 1].", call. = FALSE)
  }

  ok <- is.finite(x) & is.finite(w) & (w >= 0)
  x <- x[ok]
  w <- w[ok]

  if (length(x) == 0) return(rep(NA_real_, length(probs)))

  sw <- sum(w)
  if (!is.finite(sw) || sw <= .Machine$double.eps) {
    return(rep(NA_real_, length(probs)))
  }

  w <- w / sw
  o <- order(x)
  x <- x[o]
  w <- w[o]

  cw <- cumsum(w)

  sapply(probs, function(p) {
    x[which(cw >= p)[1]]
  })
}


# ------------------------------------------------------------------------------
# Effective ensemble size
# ------------------------------------------------------------------------------
# Neff = 1 / sum(w^2) after normalizing w to sum to 1.
# Interpretation:
# - Neff = N when weights equal
# - Neff decreases as weights concentrate
effective_n <- function(w) {
  if (!is.numeric(w)) stop("w must be numeric.", call. = FALSE)
  w <- w[is.finite(w) & w > 0]

  if (length(w) == 0) return(NA_real_)

  w <- w / sum(w)
  1 / sum(w^2)
}


# ------------------------------------------------------------------------------
# Weighted exceedance probability
# ------------------------------------------------------------------------------
# direction:
# - "gt": P(X > threshold)
# - "lt": P(X < threshold)
prob_exceed <- function(x, w, threshold, direction = c("gt", "lt")) {
  direction <- match.arg(direction)

  if (!is.numeric(x) || !is.numeric(w)) stop("x and w must be numeric.", call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1) stop("threshold must be a single numeric value.", call. = FALSE)

  ok <- is.finite(x) & is.finite(w) & (w >= 0)
  x <- x[ok]
  w <- w[ok]

  if (length(x) == 0) return(NA_real_)

  sw <- sum(w)
  if (!is.finite(sw) || sw <= .Machine$double.eps) return(NA_real_)

  w <- w / sw

  if (direction == "gt") {
    sum(w[x > threshold])
  } else {
    sum(w[x < threshold])
  }
}


# ------------------------------------------------------------------------------
# Standard reporting table per scenario
# ------------------------------------------------------------------------------
# thresholds (optional) format:
# thresholds <- list(
#   loss_gt_20 = list(value = 20, direction = "gt"),
#   deficit_lt_10 = list(value = -10, direction = "lt")
# )
report_weighted_impacts_from_grid <- function(
    grid_probs,
    stress_test_grid,
    value_col = "value",
    scenario_col = "scenario",
    probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
    thresholds = NULL) {

  # join probabilities to impacts
  df <- merge(
    grid_probs,
    stress_test_grid,
    by = c("tavg", "prcp"),
    all.x = TRUE,
    sort = FALSE
  )

  if (any(is.na(df[[value_col]]))) {
    stop("Some stress-test grid points have no impact value.", call. = FALSE)
  }

  scenarios <- sort(unique(df[[scenario_col]]))
  out <- vector("list", length(scenarios))
  names(out) <- scenarios

  for (sc in scenarios) {
    d <- df[df[[scenario_col]] == sc, , drop = FALSE]

    x <- d[[value_col]]
    w <- d$prob

    # normalize defensively
    w <- w / sum(w)

    q <- weighted_quantile(x, w, probs)
    names(q) <- paste0("P", probs * 100)

    row <- data.frame(
      scenario = sc,
      t(q),
      eff_n = effective_n(w),
      stringsAsFactors = FALSE
    )

    if (!is.null(thresholds)) {
      for (nm in names(thresholds)) {
        thr <- thresholds[[nm]]$value
        dir <- thresholds[[nm]]$direction
        row[[paste0("P_", nm)]] <- prob_exceed(x, w, thr, dir)
      }
    }

    out[[sc]] <- row
  }

  do.call(rbind, out)
}




# ------------------------------------------------------------------------------
# 1) Simple weights calculator: equal weight per institution (within each scenario)
# ------------------------------------------------------------------------------

#' Calculate simple institution-based weights (per scenario)
#'
#' @description
#' Within each scenario (e.g., SSP), assigns equal total weight to each institution,
#' then splits that institution weight equally among models from that institution
#' present in that scenario. If multiple rows exist for the same model within a
#' scenario (e.g., ensemble members), the model's weight is split equally across rows.
#'
#' @param data data.frame containing at least scenario, model, institution.
#' @param scenario_col Character. Scenario/SSP column name. Default "scenario".
#' @param institution_col Character. Institution column name. Default "institution".
#' @param model_col Character. Model column name. Default "model".
#' @param weight_col Character. Output weight column name. Default "w_inst".
#'
#' @return `data` with an added column `weight_col` summing to 1 within each scenario.
#' @export
compute_weights_institutions <- function(
    data,
    scenario_col = "scenario",
    institution_col = "institution",
    model_col = "model",
    weight_col = "w_inst") {

  required <- c(scenario_col, institution_col, model_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  df <- data
  df[[weight_col]] <- NA_real_

  scenarios <- sort(unique(as.character(df[[scenario_col]])))

  for (sc in scenarios) {
    idx_sc <- which(df[[scenario_col]] == sc)
    df_sc <- df[idx_sc, , drop = FALSE]

    # model -> institution (first non-NA)
    inst_by_model <- tapply(df_sc[[institution_col]], df_sc[[model_col]], function(z) z[which(!is.na(z))[1]])
    models_sc <- names(inst_by_model)
    insts_present <- unique(unname(inst_by_model))
    insts_present <- insts_present[!is.na(insts_present)]
    n_inst <- length(insts_present)

    if (n_inst == 0) {
      # fallback: equal per row
      df[[weight_col]][idx_sc] <- rep(1 / length(idx_sc), length(idx_sc))
      next
    }

    w_inst_total <- 1 / n_inst

    # model weights
    w_model <- numeric(length(models_sc))
    names(w_model) <- models_sc

    for (m in models_sc) {
      inst_m <- unname(inst_by_model[m])
      n_models_in_inst <- sum(unname(inst_by_model) == inst_m)
      w_model[m] <- w_inst_total / n_models_in_inst
    }

    # split model weights across rows
    for (m in models_sc) {
      idx_m <- idx_sc[df_sc[[model_col]] == m]
      df[[weight_col]][idx_m] <- w_model[m] / length(idx_m)
    }

    # renormalize (guard)
    s <- sum(df[[weight_col]][idx_sc], na.rm = TRUE)
    if (is.finite(s) && s > .Machine$double.eps) {
      df[[weight_col]][idx_sc] <- df[[weight_col]][idx_sc] / s
    }
  }

  df
}





################################################################################



#' Compute genealogy-based model weights per SSP based on Kuma et al. 2023.
#'
#' @description
#' Computes relative weights for models within each SSP scenario using the
#' genealogy information from Kuma et al. (2023) supplemental Table S1.
#'
#' Default method ("family") gives each genealogy family equal total weight
#' within each SSP, then splits that family weight equally among models from
#' that family present in that SSP. If multiple rows per model exist (e.g.,
#' different ensemble members/variants), the model weight is split equally
#' across those rows.
#'
#' @param gcm_data data.frame with at least model and scenario columns.
#' @param kuma_table data.frame read from Kuma et al. CSV (Table S1).
#' @param model_col Name of model column in `gcm_data`. Default "model".
#' @param scenario_col Name of SSP/scenario column in `gcm_data`. Default "scenario".
#' @param cmip_phase Which name column to use for mapping. One of "CMIP6", "CMIP5", "CMIP3".
#'   Default "CMIP6".
#' @param method Weighting method. One of:
#'   - "family" (default): Equal weight per family
#'   - "family_sqrt": Weight proportional to 1/sqrt(n_models_in_family)
#'   - "independence": Ignore genealogy (equal weight per model)
#' @param keep_original_model Logical. If TRUE, keep original model string and add a cleaned model field.
#' @param verbose Logical. If TRUE, prints a short mapping/coverage summary.
#'
#' @return
#' `gcm_data` with added columns:
#' - model_clean: cleaned CMIP model name used for mapping
#' - model_family: genealogy family (from Kuma table; fallback to model_clean if unmatched)
#' - w_genealogy: per-row weight within scenario_col; sums to 1 per scenario
#'
#' @export
compute_weights_genealogy <- function(
    gcm_data,
    kuma_table,
    model_col = "model",
    scenario_col = "scenario",
    cmip_phase = c("CMIP6", "CMIP5", "CMIP3"),
    method = c("family", "family_sqrt", "independence"),
    keep_original_model = TRUE,
    verbose = TRUE) {

  cmip_phase <- match.arg(cmip_phase)
  method <- match.arg(method)

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  .stop_if_missing <- function(df, cols, df_name) {
    miss <- setdiff(cols, names(df))
    if (length(miss) > 0) {
      stop(df_name, " is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
  }

  .clean_model_name <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    # handle strings like "IPSL/IPSL-CM6A-LR" (keep last token)
    x <- sub("^.*\\/", "", x)
    x <- trimws(x)
    x
  }

  .split_names <- function(x) {
    # Kuma CMIP name fields are comma-separated in the CSV (as seen for IPSL)
    if (is.na(x) || !nzchar(x)) return(character(0))
    out <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
    out <- trimws(out)
    out[nzchar(out)]
  }

  .build_mapping <- function(kuma_table, cmip_phase) {
    name_col <- switch(
      cmip_phase,
      CMIP6 = "CMIP6 names",
      CMIP5 = "CMIP5 names",
      CMIP3 = "CMIP3 names"
    )

    .stop_if_missing(kuma_table, c("Model", "Family", name_col), "kuma_table")

    # Build two-column mapping: cmip_name -> family
    cmip_name <- character(0)
    family <- character(0)

    for (i in seq_len(nrow(kuma_table))) {
      nm <- .split_names(kuma_table[[name_col]][i])
      if (length(nm) == 0) next
      cmip_name <- c(cmip_name, nm)
      family <- c(family, rep(as.character(kuma_table$Family[i]), length(nm)))
    }

    # deduplicate; if duplicates occur, keep first (should be rare)
    map <- data.frame(cmip_name = cmip_name, model_family = family, stringsAsFactors = FALSE)
    map <- map[!duplicated(map$cmip_name), , drop = FALSE]
    map
  }

  # Family-equal weights within scenario:
  # 1) for each scenario, compute weight per family = 1 / n_families_present
  # 2) split equally among models in that family in that scenario
  .compute_family_weights_per_scenario <- function(df, scenario_col) {
    df$w_genealogy <- NA_real_

    scenarios <- sort(unique(as.character(df[[scenario_col]])))

    for (sc in scenarios) {
      idx_sc <- which(df[[scenario_col]] == sc)
      df_sc <- df[idx_sc, , drop = FALSE]

      # unique models in this scenario
      models_sc <- unique(df_sc$model_clean)

      # family per model (one value per model_clean assumed)
      fam_by_model <- tapply(df_sc$model_family, df_sc$model_clean, function(z) z[1])

      fams_present <- unique(unname(fam_by_model[models_sc]))
      fams_present <- fams_present[!is.na(fams_present)]
      n_fam <- length(fams_present)

      if (n_fam == 0) {
        # Improved fallback: weight by unique models, not rows
        w_model <- 1 / length(models_sc)
        for (m in models_sc) {
          idx_m <- idx_sc[df_sc$model_clean == m]
          df$w_genealogy[idx_m] <- w_model / length(idx_m)
        }
        next
      }

      w_family <- 1 / n_fam

      # compute model weights (per unique model)
      w_model <- numeric(length(models_sc))
      names(w_model) <- models_sc

      for (m in models_sc) {
        fam_m <- unname(fam_by_model[m])
        n_models_in_fam <- sum(unname(fam_by_model[models_sc]) == fam_m)
        w_model[m] <- w_family / n_models_in_fam
      }

      # split model weight across multiple rows per model in this scenario (members/variants)
      for (m in models_sc) {
        idx_m <- idx_sc[df_sc$model_clean == m]
        n_rows_m <- length(idx_m)
        df$w_genealogy[idx_m] <- w_model[m] / n_rows_m
      }

      # numerical guard: renormalize within scenario
      s <- sum(df$w_genealogy[idx_sc], na.rm = TRUE)
      if (is.finite(s) && s > .Machine$double.eps) {
        df$w_genealogy[idx_sc] <- df$w_genealogy[idx_sc] / s
      }
    }

    df
  }

  # Family-sqrt weights within scenario:
  # Weight proportional to 1/sqrt(n_models_in_family)
  .compute_family_sqrt_weights_per_scenario <- function(df, scenario_col) {
    df$w_genealogy <- NA_real_

    scenarios <- sort(unique(as.character(df[[scenario_col]])))

    for (sc in scenarios) {
      idx_sc <- which(df[[scenario_col]] == sc)
      df_sc <- df[idx_sc, , drop = FALSE]

      # unique models in this scenario
      models_sc <- unique(df_sc$model_clean)

      # family per model (one value per model_clean assumed)
      fam_by_model <- tapply(df_sc$model_family, df_sc$model_clean, function(z) z[1])

      fams_present <- unique(unname(fam_by_model[models_sc]))
      fams_present <- fams_present[!is.na(fams_present)]
      n_fam <- length(fams_present)

      if (n_fam == 0) {
        # Fallback: equal model weights
        w_model <- 1 / length(models_sc)
        for (m in models_sc) {
          idx_m <- idx_sc[df_sc$model_clean == m]
          df$w_genealogy[idx_m] <- w_model / length(idx_m)
        }
        next
      }

      # compute family weights: proportional to 1/sqrt(n_models_in_family)
      family_weights <- numeric(length(fams_present))
      names(family_weights) <- fams_present

      for (fam in fams_present) {
        n_models_in_fam <- sum(unname(fam_by_model[models_sc]) == fam)
        family_weights[fam] <- 1 / sqrt(n_models_in_fam)
      }

      # normalize family weights to sum to 1
      family_weights <- family_weights / sum(family_weights)

      # compute model weights (per unique model)
      w_model <- numeric(length(models_sc))
      names(w_model) <- models_sc

      for (m in models_sc) {
        fam_m <- unname(fam_by_model[m])
        n_models_in_fam <- sum(unname(fam_by_model[models_sc]) == fam_m)
        w_model[m] <- family_weights[fam_m] / n_models_in_fam
      }

      # split model weight across multiple rows per model in this scenario (members/variants)
      for (m in models_sc) {
        idx_m <- idx_sc[df_sc$model_clean == m]
        n_rows_m <- length(idx_m)
        df$w_genealogy[idx_m] <- w_model[m] / n_rows_m
      }

      # numerical guard: renormalize within scenario
      s <- sum(df$w_genealogy[idx_sc], na.rm = TRUE)
      if (is.finite(s) && s > .Machine$double.eps) {
        df$w_genealogy[idx_sc] <- df$w_genealogy[idx_sc] / s
      }
    }

    df
  }

  # Independence weights within scenario:
  # Equal weight per model (ignores genealogy)
  .compute_independence_weights_per_scenario <- function(df, scenario_col) {
    df$w_genealogy <- NA_real_

    scenarios <- sort(unique(as.character(df[[scenario_col]])))

    for (sc in scenarios) {
      idx_sc <- which(df[[scenario_col]] == sc)
      df_sc <- df[idx_sc, , drop = FALSE]

      # unique models in this scenario
      models_sc <- unique(df_sc$model_clean)
      w_model <- 1 / length(models_sc)

      # split model weight across multiple rows per model in this scenario
      for (m in models_sc) {
        idx_m <- idx_sc[df_sc$model_clean == m]
        df$w_genealogy[idx_m] <- w_model / length(idx_m)
      }

      # numerical guard: renormalize within scenario
      s <- sum(df$w_genealogy[idx_sc], na.rm = TRUE)
      if (is.finite(s) && s > .Machine$double.eps) {
        df$w_genealogy[idx_sc] <- df$w_genealogy[idx_sc] / s
      }
    }

    df
  }

  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------

  .stop_if_missing(gcm_data, c(model_col, scenario_col), "gcm_data")

  # ---------------------------------------------------------------------------
  # Build mapping and attach families (single pass - Suggestion #1)
  # ---------------------------------------------------------------------------

  map <- .build_mapping(kuma_table, cmip_phase = cmip_phase)

  # Work on a copy with row ID for stable ordering
  df <- gcm_data
  df$.row_id <- seq_len(nrow(df))

  if (keep_original_model && model_col != "model_raw") {
    df$model_raw <- df[[model_col]]
  }

  df$model_clean <- .clean_model_name(df[[model_col]])

  # attach model_family
  df <- merge(
    df,
    map,
    by.x = "model_clean",
    by.y = "cmip_name",
    all.x = TRUE,
    sort = FALSE
  )

  # fallback for unmatched: treat as its own family
  unmatched <- is.na(df$model_family) | !nzchar(df$model_family)
  if (any(unmatched)) {
    df$model_family[unmatched] <- df$model_clean[unmatched]
  }

  # restore original row order
  df <- df[order(df$.row_id), , drop = FALSE]
  df$.row_id <- NULL

  # ---------------------------------------------------------------------------
  # Compute weights (Suggestion #4 - multiple methods)
  # ---------------------------------------------------------------------------

  if (method == "family") {
    df <- .compute_family_weights_per_scenario(df, scenario_col = scenario_col)
  } else if (method == "family_sqrt") {
    df <- .compute_family_sqrt_weights_per_scenario(df, scenario_col = scenario_col)
  } else if (method == "independence") {
    df <- .compute_independence_weights_per_scenario(df, scenario_col = scenario_col)
  }

  # ---------------------------------------------------------------------------
  # Messaging (Suggestion #3 - enhanced diagnostics)
  # ---------------------------------------------------------------------------

  if (verbose) {
    n_total_models <- length(unique(df$model_clean))
    n_mapped <- sum(unique(df$model_clean) %in% map$cmip_name)
    n_unmapped <- n_total_models - n_mapped
    message(sprintf(
      "[GENEALOGY] CMIP phase=%s | method=%s | unique models=%d | mapped=%d | unmapped=%d (fallback family=model)",
      cmip_phase, method, n_total_models, n_mapped, n_unmapped
    ))

    # Per-scenario diagnostics
    scenarios <- sort(unique(as.character(df[[scenario_col]])))
    for (sc in scenarios) {
      df_sc <- df[df[[scenario_col]] == sc, , drop = FALSE]
      n_fam_sc <- length(unique(df_sc$model_family))
      n_mod_sc <- length(unique(df_sc$model_clean))
      message(sprintf("  %s: %d families | %d models", sc, n_fam_sc, n_mod_sc))
    }
  }

  df
}
