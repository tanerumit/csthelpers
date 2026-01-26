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
