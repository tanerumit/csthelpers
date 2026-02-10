# =============================================================================
# Performance–Independence Weighting (PIW) for Climate Model Ensembles
# Drop-in replacement script (patched methodological + technical issues)
#
# Main functions:
#   - compute_gcm_weights_bma()
#   - compute_weighted_ensemble_stats()
#
# Notes:
# - Base-R only in this script (assumes external helpers exist in your package):
#     compute_gcm_weights_by_genealogy()
#     compute_gcm_weights_by_institution()
# =============================================================================

# -----------------------------------------------------------------------------
# Helpers (internal)
# -----------------------------------------------------------------------------

#' Null-coalescing operator (internal)
#' @description Returns `a` if not NULL, otherwise `b`.
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Robust scale estimator (internal)
#' @description Returns a robust scale estimate (sd or mad) with safe fallbacks.
#' @keywords internal
.piw_scale <- function(x, robust_scale = c("sd", "mad"), min_scale = 1e-8) {
  robust_scale <- match.arg(robust_scale)
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(1)

  s <- switch(
    robust_scale,
    "sd"  = stats::sd(x),
    "mad" = stats::mad(x, constant = 1) # constant=1 => MAD in same units
  )
  if (!is.finite(s) || s < min_scale) return(1)
  s
}

#' Map skill metric to a score where higher is better (internal)
#' @description Converts a skill metric to a "higher is better" score for weighting.
#' @keywords internal
.skill_to_score <- function(skill, metric) {
  if (metric %in% c("rmse", "bias")) return(-abs(skill))
  if (metric %in% c("correlation", "nse")) return(skill)
  stop("Unknown 'metric': ", metric, call. = FALSE)
}

#' Safe column resolver (internal)
#' @description Returns `df[[primary]]` if present else `df[[fallback]]`.
#' @keywords internal
.piw_col_or <- function(df, primary, fallback) {
  if (!is.data.frame(df)) stop(".piw_col_or(): df must be a data.frame.", call. = FALSE)
  if (primary %in% names(df)) return(df[[primary]])
  if (fallback %in% names(df)) return(df[[fallback]])
  stop("Missing expected columns: '", primary, "' and fallback '", fallback, "'.", call. = FALSE)
}

#' Entropy of a weight vector (internal)
#' @description Computes Shannon entropy for a nonnegative weight vector.
#' @keywords internal
.piw_entropy <- function(w) {
  w <- w[is.finite(w)]
  w <- w[w > 0]
  if (length(w) == 0) return(NA_real_)
  -sum(w * log(w))
}

#' Gini coefficient of a weight vector (internal)
#' @description Computes a Gini coefficient for a nonnegative weight vector.
#' @keywords internal
.piw_gini <- function(w) {
  w <- w[is.finite(w)]
  if (length(w) == 0) return(NA_real_)
  w <- pmax(w, 0)
  s <- sum(w)
  if (s < .Machine$double.eps) return(NA_real_)
  w <- w / s

  w_sorted <- sort(w)
  n <- length(w_sorted)
  1 - 2 * sum((n:1) * w_sorted) / n + 1 / n
}

# -----------------------------------------------------------------------------
# Helpers (semi-public building blocks; not exported unless your NAMESPACE exports)
# -----------------------------------------------------------------------------

#' Convert model skill to performance weights (internal helper)
#' @description Converts model skill values into normalized decision weights.
#' @keywords internal
skill_to_weight <- function(
    skill,
    skill_metric = c("rmse", "correlation", "bias", "nse"),
    transform = c("gaussian", "softmax", "rank", "linear"),
    temperature = NULL,
    robust_scale = c("sd", "mad"),
    epsilon = 1e-3
) {
  skill_metric <- match.arg(skill_metric)
  transform <- match.arg(transform)
  robust_scale <- match.arg(robust_scale)

  score <- .skill_to_score(skill, skill_metric)

  finite_idx <- which(is.finite(score))
  out <- rep(NA_real_, length(score))
  if (length(finite_idx) == 0) return(out)

  score_f <- score[finite_idx]
  if (length(unique(score_f)) <= 1) {
    out[finite_idx] <- rep(1 / length(finite_idx), length(finite_idx))
    return(out)
  }

  if (is.null(temperature)) {
    temperature <- .piw_scale(score_f, robust_scale = robust_scale)
  }
  if (!is.numeric(temperature) || length(temperature) != 1 || !is.finite(temperature) || temperature <= 0) {
    stop("'temperature' must be a positive finite numeric scalar (or NULL).", call. = FALSE)
  }

  w_raw <- switch(
    transform,
    "softmax" = {
      s0 <- score_f - max(score_f)
      exp(s0 / temperature)
    },
    "gaussian" = {
      d <- (max(score_f) - score_f) / temperature
      exp(-0.5 * d^2)
    },
    "rank" = {
      rank(score_f, ties.method = "average")
    },
    "linear" = {
      s_min <- min(score_f)
      s_max <- max(score_f)
      rng <- s_max - s_min
      if (!is.finite(rng) || rng < .Machine$double.eps) {
        rep(1, length(score_f))
      } else {
        (score_f - s_min) / rng + epsilon
      }
    }
  )

  s <- sum(w_raw)
  if (!is.finite(s) || s < .Machine$double.eps) {
    out[finite_idx] <- rep(1 / length(finite_idx), length(finite_idx))
    return(out)
  }

  out[finite_idx] <- w_raw / s
  out
}

#' Compute a single skill metric (internal helper)
#' @description Computes a scalar skill metric for aligned model vs observed values.
#' @keywords internal
compute_skill_metric <- function(
    gcm_values,
    obs_values,
    skill_metric = c("rmse", "correlation", "bias", "nse")
) {
  skill_metric <- match.arg(skill_metric)

  gcm_values <- as.numeric(gcm_values)
  obs_values <- as.numeric(obs_values)

  if (length(gcm_values) != length(obs_values)) {
    stop("compute_skill_metric(): gcm_values and obs_values must have identical length after alignment.", call. = FALSE)
  }

  ok <- is.finite(gcm_values) & is.finite(obs_values)
  if (!any(ok)) return(NA_real_)

  g <- gcm_values[ok]
  o <- obs_values[ok]

  switch(
    skill_metric,
    "rmse" = sqrt(mean((g - o)^2)),
    "bias" = mean(g - o),
    "correlation" = {
      if (stats::sd(g) < .Machine$double.eps || stats::sd(o) < .Machine$double.eps) return(0)
      stats::cor(g, o)
    },
    "nse" = {
      ss_res <- sum((o - g)^2)
      ss_tot <- sum((o - mean(o))^2)
      if (!is.finite(ss_tot) || ss_tot < .Machine$double.eps) return(0)
      1 - ss_res / ss_tot
    }
  )
}

#' Combine structural and performance weights (internal helper)
#' @description Combines two weight vectors using a log-space geometric blend.
#' @keywords internal
combine_weights_piw <- function(w_structure, w_skill, skill_weight, min_weight = 1e-3) {
  if (length(w_structure) != length(w_skill)) stop("Weight vectors must have equal length.", call. = FALSE)
  if (!is.numeric(skill_weight) || length(skill_weight) != 1 || !is.finite(skill_weight) || skill_weight < 0 || skill_weight > 1) {
    stop("'skill_weight' must be a finite numeric scalar in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(min_weight) || length(min_weight) != 1 || !is.finite(min_weight) || min_weight <= 0 || min_weight >= 1) {
    stop("'min_weight' must be a finite numeric scalar in (0, 1).", call. = FALSE)
  }

  w_s <- pmax(w_structure, min_weight)
  w_k <- pmax(w_skill, min_weight)

  lp <- (1 - skill_weight) * log(w_s) + skill_weight * log(w_k)
  lp <- lp - max(lp)
  w <- exp(lp)
  w / sum(w)
}

# =============================================================================
# Main functions
# =============================================================================

#' Compute PIW (Performance–Independence Weighting) ensemble weights
#'
#' @description
#' Computes decision weights for climate model ensembles by combining a
#' structural (independence) component with a performance (skill) component.
#'
#' @details
#' The function estimates model-level skill (aggregated across scenarios by
#' default) and applies that performance weighting across scenarios, then blends
#' with structural weights per scenario via a log-space geometric mixture
#' controlled by `skill_weight`.
#'
#' Structural priors require external package helpers for genealogy/institution
#' priors; uniform structural weights are computed internally.
#'
#' @param gcm_data data.frame containing at least `model_col` and `scenario_col`.
#' @param skill_data Optional data.frame with per-model (or per-model-scenario)
#'   skill in `skill_col`.
#' @param obs_data Optional data.frame of observations used to compute skill
#'   when `skill_data` is NULL.
#' @param kuma_table Required for genealogy priors; passed to
#'   `compute_gcm_weights_by_genealogy()`.
#' @param model_col Character scalar column name identifying model.
#' @param scenario_col Character scalar column name identifying scenario.
#' @param skill_col Character scalar skill column name in `skill_data`.
#' @param target_col Character scalar target column name in `gcm_data`/`obs_data`
#'   required when computing skill from observations.
#' @param prior_method Structural weighting method: `"genealogy"`,
#'   `"genealogy_sqrt"`, `"institution"`, or `"uniform"`.
#' @param institution_col Character scalar institution column in `gcm_data`
#'   required when `prior_method = "institution"`.
#' @param skill_metric Skill metric used for performance weighting: `"rmse"`,
#'   `"correlation"`, `"bias"`, or `"nse"`.
#' @param likelihood_transform Preference transform mapping skill to weights:
#'   `"gaussian"`, `"softmax"`, `"rank"`, or `"linear"`.
#' @param skill_weight Numeric scalar in [0,1] controlling the blend between
#'   structure-only (0) and skill-only (1).
#' @param min_weight Numeric scalar in (0,1) used to regularize near-zero weights.
#' @param verbose Logical; if TRUE prints progress and summary diagnostics.
#' @param temperature Optional numeric scalar > 0 controlling transform sharpness;
#'   if NULL it is estimated from skill dispersion using `robust_scale`.
#' @param robust_scale Robust scale estimator used when `temperature` is NULL:
#'   `"sd"` or `"mad"`.
#' @param include_missing_skill Logical; if FALSE drops models without finite
#'   skill values, else retains them with minimal regularized influence.
#' @param skill_weight_grid Optional numeric vector of values in [0,1] used to
#'   compute sensitivity diagnostics; NULL disables sensitivity output.
#'
#' @return
#' A data.frame with columns:
#' \itemize{
#'   \item model, scenario
#'   \item w_prior (structural), w_likelihood (skill), w_bma (PIW combined)
#'   \item skill_value, w_structure_only, w_skill_only
#' }
#' Additional diagnostics are attached as attributes:
#' \itemize{
#'   \item attr(result, "piw_diagnostics")
#'   \item attr(result, "piw_sensitivity")
#' }
#'
#' @examples
#' \dontrun{
#' w <- compute_gcm_weights_bma(
#'   gcm_data = gcm,
#'   skill_data = skill,
#'   prior_method = "institution",
#'   skill_metric = "rmse",
#'   likelihood_transform = "gaussian",
#'   skill_weight = 0.5
#' )
#' attr(w, "piw_diagnostics")
#' }
#'
#' @export
compute_gcm_weights_bma <- function(
    gcm_data,
    skill_data = NULL,
    obs_data = NULL,
    kuma_table = NULL,
    model_col = "model",
    scenario_col = "scenario",
    skill_col = "skill",
    target_col = NULL,
    prior_method = c("genealogy", "genealogy_sqrt", "institution", "uniform"),
    institution_col = "institution",
    skill_metric = c("rmse", "correlation", "bias", "nse"),
    likelihood_transform = c("gaussian", "softmax", "rank", "linear"),
    skill_weight = 0.5,
    min_weight = 0.001,
    verbose = TRUE,
    temperature = NULL,
    robust_scale = c("sd", "mad"),
    include_missing_skill = FALSE,
    skill_weight_grid = seq(0, 1, by = 0.25)
) {
  prior_method <- match.arg(prior_method)
  skill_metric <- match.arg(skill_metric)
  likelihood_transform <- match.arg(likelihood_transform)
  robust_scale <- match.arg(robust_scale)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(gcm_data)) stop("'gcm_data' must be a data.frame.", call. = FALSE)
  if (!model_col %in% names(gcm_data)) stop("Column '", model_col, "' not found in gcm_data.", call. = FALSE)
  if (!scenario_col %in% names(gcm_data)) stop("Column '", scenario_col, "' not found in gcm_data.", call. = FALSE)

  if (!is.numeric(skill_weight) || length(skill_weight) != 1 || !is.finite(skill_weight) || skill_weight < 0 || skill_weight > 1) {
    stop("'skill_weight' must be a numeric scalar in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(min_weight) || length(min_weight) != 1 || !is.finite(min_weight) || min_weight <= 0 || min_weight >= 1) {
    stop("'min_weight' must be a numeric scalar in (0, 1).", call. = FALSE)
  }
  if (!is.logical(include_missing_skill) || length(include_missing_skill) != 1 || is.na(include_missing_skill)) {
    stop("'include_missing_skill' must be TRUE/FALSE.", call. = FALSE)
  }

  if (is.null(skill_data) && is.null(obs_data)) {
    stop("Either 'skill_data' or 'obs_data' must be provided.", call. = FALSE)
  }
  if (prior_method %in% c("genealogy", "genealogy_sqrt") && is.null(kuma_table)) {
    stop("'kuma_table' is required for prior_method = '", prior_method, "'.", call. = FALSE)
  }
  if (prior_method == "institution" && !institution_col %in% names(gcm_data)) {
    stop("Column '", institution_col, "' not found (required for institution prior).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 1: Model-level skill table (decoupled from scenarios)
  # ---------------------------------------------------------------------------
  if (is.null(skill_data)) {
    if (isTRUE(verbose)) message("[PIW] Computing skill scores from observations (requires pre-aligned data).")

    if (is.null(target_col) || !target_col %in% names(gcm_data)) {
      stop("'target_col' must be specified and present in gcm_data to compute skill.", call. = FALSE)
    }
    if (!is.data.frame(obs_data) || !target_col %in% names(obs_data)) {
      stop("'obs_data' must be a data.frame containing 'target_col'.", call. = FALSE)
    }

    models <- unique(gcm_data[[model_col]])
    skill_model <- data.frame(
      model = models,
      skill_value = NA_real_,
      stringsAsFactors = FALSE
    )

    obs_vals <- obs_data[[target_col]]

    for (i in seq_along(models)) {
      m <- models[i]
      mod_vals <- gcm_data[gcm_data[[model_col]] == m, target_col, drop = TRUE]

      n <- min(length(mod_vals), length(obs_vals))
      if (n <= 0) next

      sv <- compute_skill_metric(
        gcm_values = mod_vals[seq_len(n)],
        obs_values = obs_vals[seq_len(n)],
        skill_metric = skill_metric
      )
      skill_model$skill_value[i] <- sv
    }
  } else {
    if (!is.data.frame(skill_data)) stop("'skill_data' must be a data.frame.", call. = FALSE)
    if (!model_col %in% names(skill_data)) stop("Column '", model_col, "' not found in skill_data.", call. = FALSE)
    if (!skill_col %in% names(skill_data)) stop("Column '", skill_col, "' not found in skill_data.", call. = FALSE)

    if (isTRUE(verbose) && scenario_col %in% names(skill_data)) {
      message("[PIW] skill_data contains scenarios; PIW aggregates skill at model level and applies to all scenarios.")
    }

    models <- unique(skill_data[[model_col]])
    skill_model <- data.frame(
      model = models,
      skill_value = NA_real_,
      stringsAsFactors = FALSE
    )

    for (i in seq_along(models)) {
      m <- models[i]
      vals <- skill_data[[skill_col]][skill_data[[model_col]] == m]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) next
      skill_model$skill_value[i] <- mean(vals)
    }
  }

  if (!include_missing_skill) {
    keep <- is.finite(skill_model$skill_value)
    if (isTRUE(verbose)) message("[PIW] Dropping models with missing skill: ", sum(!keep), " dropped.")
    skill_model <- skill_model[keep, , drop = FALSE]
  } else {
    if (isTRUE(verbose)) message("[PIW] include_missing_skill=TRUE: missing skill will be regularized via min_weight.")
  }

  # ---------------------------------------------------------------------------
  # Step 2: Structural weights (per scenario)
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) message("[PIW] Computing structural weights (method: ", prior_method, ").")

  if (prior_method == "genealogy") {
    prior_df <- compute_gcm_weights_by_genealogy(
      gcm_data = gcm_data,
      kuma_table = kuma_table,
      model_col = model_col,
      scenario_col = scenario_col,
      method = "family",
      verbose = FALSE
    )
    prior_df$w_prior <- prior_df$w_genealogy
  } else if (prior_method == "genealogy_sqrt") {
    prior_df <- compute_gcm_weights_by_genealogy(
      gcm_data = gcm_data,
      kuma_table = kuma_table,
      model_col = model_col,
      scenario_col = scenario_col,
      method = "family_sqrt",
      verbose = FALSE
    )
    prior_df$w_prior <- prior_df$w_genealogy
  } else if (prior_method == "institution") {
    prior_df <- compute_gcm_weights_by_institution(
      data = gcm_data,
      scenario_col = scenario_col,
      institution_col = institution_col,
      model_col = model_col,
      weight_col = "w_prior"
    )
  } else {
    ms <- unique(gcm_data[, c(model_col, scenario_col), drop = FALSE])
    names(ms) <- c("model", "scenario")

    ms$w_prior <- NA_real_
    scenarios <- unique(ms$scenario)

    for (sc in scenarios) {
      idx <- which(ms$scenario == sc)
      models_sc <- unique(ms$model[idx])
      if (length(models_sc) == 0) next
      ms$w_prior[idx] <- 1 / length(models_sc)
    }
    prior_df <- ms
  }

  # Aggregate structural weights to model-scenario safely.
  prior_agg <- stats::aggregate(
    prior_df["w_prior"],
    by = list(
      model = .piw_col_or(prior_df, "model", model_col),
      scenario = .piw_col_or(prior_df, "scenario", scenario_col)
    ),
    FUN = sum,
    na.rm = TRUE
  )
  names(prior_agg)[1:2] <- c("model", "scenario")

  # ---------------------------------------------------------------------------
  # Step 3: Performance weights (model-level)
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    msg_t <- if (is.null(temperature)) "auto" else format(temperature)
    message(
      "[PIW] Computing performance weights (metric: ", skill_metric,
      ", transform: ", likelihood_transform,
      ", temperature: ", msg_t, ", robust_scale: ", robust_scale, ")."
    )
  }

  perf_w <- skill_to_weight(
    skill = skill_model$skill_value,
    skill_metric = skill_metric,
    transform = likelihood_transform,
    temperature = temperature,
    robust_scale = robust_scale
  )

  perf_df <- data.frame(
    model = skill_model$model,
    skill_value = skill_model$skill_value,
    w_likelihood = perf_w,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Step 4: Merge + combine to PIW weights per scenario
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) message("[PIW] Combining structural + performance (skill_weight: ", skill_weight, ").")

  result_df <- merge(prior_agg, perf_df, by = "model", all.x = TRUE, all.y = FALSE)

  if (include_missing_skill) {
    result_df$w_likelihood[!is.finite(result_df$w_likelihood)] <- min_weight
    result_df$skill_value[!is.finite(result_df$skill_value)] <- NA_real_
  } else {
    keep <- is.finite(result_df$w_likelihood)
    if (isTRUE(verbose)) message("[PIW] Dropping model-scenario rows lacking performance weights: ", sum(!keep), " dropped.")
    result_df <- result_df[keep, , drop = FALSE]
  }

  result_df$w_prior[!is.finite(result_df$w_prior)] <- min_weight

  result_df$w_bma <- NA_real_ # retained output name for backward compatibility
  result_df$w_structure_only <- NA_real_
  result_df$w_skill_only <- NA_real_

  scenarios <- unique(result_df$scenario)
  diag_list <- list()

  for (sc in scenarios) {
    idx <- which(result_df$scenario == sc)
    if (length(idx) == 0) next

    w_s <- result_df$w_prior[idx]
    w_k <- result_df$w_likelihood[idx]

    w_s0 <- pmax(w_s, min_weight)
    w_s0 <- w_s0 / sum(w_s0)

    w_k0 <- pmax(w_k, min_weight)
    w_k0 <- w_k0 / sum(w_k0)

    result_df$w_structure_only[idx] <- w_s0
    result_df$w_skill_only[idx] <- w_k0

    w_piw <- combine_weights_piw(
      w_structure = w_s0,
      w_skill = w_k0,
      skill_weight = skill_weight,
      min_weight = min_weight
    )
    result_df$w_bma[idx] <- w_piw

    n_models <- length(idx)
    n_eff <- 1 / sum(w_piw^2)
    ent <- .piw_entropy(w_piw)
    gini <- .piw_gini(w_piw)

    top_local <- which.max(w_piw)
    top_model <- result_df$model[idx[top_local]]
    top_weight <- w_piw[top_local]

    diag_list[[as.character(sc)]] <- list(
      scenario = sc,
      n_models = n_models,
      n_eff = n_eff,
      entropy = ent,
      gini = gini,
      top_model = top_model,
      top_weight = top_weight
    )
  }

  # ---------------------------------------------------------------------------
  # Step 5: Sensitivity grid (diagnostic only)
  # ---------------------------------------------------------------------------
  sens_df <- NULL
  if (!is.null(skill_weight_grid)) {
    swg <- as.numeric(skill_weight_grid)
    swg <- swg[is.finite(swg) & swg >= 0 & swg <= 1]
    swg <- unique(swg)

    if (length(swg) > 0) {
      sens_rows <- list()

      for (sc in scenarios) {
        idx <- which(result_df$scenario == sc)
        if (length(idx) == 0) next

        w_s0 <- result_df$w_structure_only[idx]
        w_k0 <- result_df$w_skill_only[idx]

        for (sw in swg) {
          w_sw <- combine_weights_piw(
            w_structure = w_s0,
            w_skill = w_k0,
            skill_weight = sw,
            min_weight = min_weight
          )
          sens_rows[[length(sens_rows) + 1]] <- data.frame(
            scenario = sc,
            skill_weight = sw,
            n_models = length(w_sw),
            n_eff = 1 / sum(w_sw^2),
            entropy = .piw_entropy(w_sw),
            gini = .piw_gini(w_sw),
            stringsAsFactors = FALSE
          )
        }
      }

      sens_df <- do.call(rbind, sens_rows)
    }
  }

  attr(result_df, "piw_diagnostics") <- diag_list
  attr(result_df, "piw_sensitivity") <- sens_df

  if (isTRUE(verbose)) {
    message("\n[PIW] Weight summary by scenario:")
    for (sc in names(diag_list)) {
      d <- diag_list[[sc]]
      message(sprintf(
        "  %s: n_models=%d, N_eff=%.1f, entropy=%.2f, gini=%.2f, top=%s (%.1f%%)",
        d$scenario, d$n_models, d$n_eff, d$entropy, d$gini, d$top_model, 100 * d$top_weight
      ))
    }
    if (!is.null(sens_df)) {
      message("\n[PIW] Sensitivity grid computed. Access via: attr(result, 'piw_sensitivity').")
    }
  }

  result_df
}

# -----------------------------------------------------------------------------
# Weighted ensemble statistics
# -----------------------------------------------------------------------------

#' Compute weighted ensemble summary statistics by scenario
#'
#' @description
#' Computes weighted mean, weighted standard deviation, effective sample size,
#' and weighted quantiles per scenario using model weights.
#'
#' @details
#' Weights are merged into `gcm_data` at the model-scenario level. Non-finite
#' weights are set to zero (explicit exclusion). If the remaining weight sum is
#' degenerate, the function falls back to uniform weights over remaining records.
#'
#' Weighted quantiles are computed by CDF inversion over sorted values.
#'
#' @param gcm_data data.frame containing model projections including `target_col`,
#'   `model_col`, and `scenario_col`.
#' @param weights_df data.frame containing weights with columns `model`, `scenario`,
#'   and `weight_col`.
#' @param model_col Character scalar model column name in `gcm_data`.
#' @param scenario_col Character scalar scenario column name in `gcm_data`.
#' @param target_col Character scalar target variable column in `gcm_data`.
#' @param weight_col Character scalar weight column in `weights_df`. Default `"w_bma"`.
#' @param quantiles Numeric vector of probabilities in [0,1] for weighted quantiles.
#'
#' @return data.frame with one row per scenario and columns:
#' `scenario`, `weighted_mean`, `weighted_sd`, `n_records`, `n_eff`, and quantiles.
#'
#' @examples
#' \dontrun{
#' stats <- compute_weighted_ensemble_stats(
#'   gcm_data = gcm,
#'   weights_df = w,
#'   target_col = "flow"
#' )
#' }
#'
#' @export
compute_weighted_ensemble_stats <- function(
    gcm_data,
    weights_df,
    model_col = "model",
    scenario_col = "scenario",
    target_col,
    weight_col = "w_bma",
    quantiles = c(0.1, 0.5, 0.9)
) {
  if (!is.data.frame(gcm_data)) stop("'gcm_data' must be a data.frame.", call. = FALSE)
  if (!is.data.frame(weights_df)) stop("'weights_df' must be a data.frame.", call. = FALSE)
  if (!target_col %in% names(gcm_data)) stop("Column '", target_col, "' not found in gcm_data.", call. = FALSE)
  if (!weight_col %in% names(weights_df)) stop("Column '", weight_col, "' not found in weights_df.", call. = FALSE)
  if (!all(c(model_col, scenario_col) %in% names(gcm_data))) stop("gcm_data missing model/scenario columns.", call. = FALSE)

  merged <- merge(
    gcm_data,
    weights_df[, c("model", "scenario", weight_col)],
    by.x = c(model_col, scenario_col),
    by.y = c("model", "scenario"),
    all.x = TRUE
  )

  merged[[weight_col]][!is.finite(merged[[weight_col]])] <- 0

  scenarios <- unique(merged[[scenario_col]])
  results <- list()

  for (sc in scenarios) {
    df_sc <- merged[merged[[scenario_col]] == sc, , drop = FALSE]
    values <- df_sc[[target_col]]
    weights <- df_sc[[weight_col]]

    ok <- is.finite(values) & is.finite(weights) & weights >= 0
    values <- values[ok]
    weights <- weights[ok]

    if (length(values) == 0) next

    w_sum <- sum(weights)
    if (!is.finite(w_sum) || w_sum < .Machine$double.eps) {
      weights <- rep(1 / length(weights), length(weights))
    } else {
      weights <- weights / w_sum
    }

    w_mean <- sum(values * weights)
    w_var <- sum(weights * (values - w_mean)^2)
    w_sd <- sqrt(max(w_var, 0))

    ord <- order(values)
    v <- values[ord]
    w <- weights[ord]
    cw <- cumsum(w)

    w_quantiles <- sapply(quantiles, function(q) {
      if (!is.finite(q) || q < 0 || q > 1) return(NA_real_)

      idx <- which(cw >= q)[1]
      if (is.na(idx)) return(v[length(v)])
      if (idx == 1) return(v[1])

      w_below <- cw[idx - 1]
      w_at <- cw[idx]
      if (!is.finite(w_at - w_below) || (w_at - w_below) < .Machine$double.eps) return(v[idx])

      frac <- (q - w_below) / (w_at - w_below)
      v[idx - 1] + frac * (v[idx] - v[idx - 1])
    })
    names(w_quantiles) <- paste0("q", quantiles * 100)

    results[[as.character(sc)]] <- data.frame(
      scenario = sc,
      weighted_mean = w_mean,
      weighted_sd = w_sd,
      n_records = length(values),
      n_eff = 1 / sum(weights^2),
      as.list(w_quantiles),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}
