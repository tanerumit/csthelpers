# =============================================================================
# GCM Weighting Functions (Improved)
# =============================================================================
#
# Two functions for computing climate model ensemble weights:
# - compute_gcm_weights_by_institution(): Institution-based weighting
# - compute_gcm_weights_by_genealogy(): Genealogy-based weighting (Kuma et al.)
#
# Improvements over original:
# 1. Vectorized operations using dplyr (O(n) instead of O(n²))
# 2. Proper tidy evaluation with .data pronoun
# 3. Configurable output column names
# 4. Weight sum validation diagnostics
# 5. Defensive checks for edge cases
# 6. Safer tapply/split handling
#
# Dependencies: dplyr, rlang
# =============================================================================

#' @importFrom dplyr group_by mutate ungroup n across all_of left_join distinct
#' @importFrom rlang .data sym

# =============================================================================
# compute_gcm_weights_by_institution
# =============================================================================

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
#' @param verbose Logical. If TRUE, prints diagnostics and validates weight sums.
#'
#' @return `data` with an added column `weight_col` summing to 1 within each scenario.
#'
#' @examples
#' \dontrun{
#' gcm_data <- data.frame(
#'   model = c("ACCESS-CM2", "ACCESS-ESM1-5", "GFDL-CM4", "IPSL-CM6A-LR"),
#'   institution = c("CSIRO", "CSIRO", "NOAA-GFDL", "IPSL"),
#'   scenario = rep("SSP2-4.5", 4)
#' )
#' result <- compute_gcm_weights_by_institution(gcm_data)
#' }
#'
#' @export
compute_gcm_weights_by_institution <- function(
    data,
    scenario_col = "scenario",
    institution_col = "institution",
    model_col = "model",
    weight_col = "w_inst",
    verbose = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }

  required <- c(scenario_col, institution_col, model_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (nrow(data) == 0) {
    warning("Input data has zero rows.", call. = FALSE)
    data[[weight_col]] <- numeric(0)
    return(data)
  }

  # ---------------------------------------------------------------------------
  # Vectorized weight computation using dplyr
  # ---------------------------------------------------------------------------

  # Step 1: Build model -> institution mapping (first non-NA institution per model)
  model_inst_map <- data |>
    dplyr::select(
      dplyr::all_of(c(scenario_col, model_col, institution_col))
    ) |>
    dplyr::filter(!is.na(.data[[institution_col]])) |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(c(scenario_col, model_col)))
    ) |>
    dplyr::summarise(
      .inst = .data[[institution_col]][1],
      .groups = "drop"
    )

  # Step 2: Count models per institution per scenario
  inst_model_counts <- model_inst_map |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(scenario_col)),
      .data$.inst
    ) |>
    dplyr::summarise(
      .n_models_in_inst = dplyr::n(),
      .groups = "drop"
    )

  # Step 3: Count institutions per scenario
  inst_counts <- inst_model_counts |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenario_col))) |>
    dplyr::summarise(
      .n_inst = dplyr::n(),
      .groups = "drop"
    )

  # Step 4: Compute model weights
  model_weights <- model_inst_map |>
    dplyr::left_join(inst_counts, by = scenario_col) |>
    dplyr::left_join(inst_model_counts, by = c(scenario_col, ".inst")) |>
    dplyr::mutate(
      .w_model = (1 / .data$.n_inst) / .data$.n_models_in_inst
    ) |>
    dplyr::select(dplyr::all_of(c(scenario_col, model_col, ".w_model")))

  # Step 5: Join weights back and split across rows
  df <- data
  df$.row_id <- seq_len(nrow(df))

  df <- df |>
    dplyr::left_join(model_weights, by = c(scenario_col, model_col)) |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(c(scenario_col, model_col)))
    ) |>
    dplyr::mutate(
      !!weight_col := .data$.w_model / dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-".w_model")

  # Handle any NA weights (models with no institution mapping)
  na_weights <- is.na(df[[weight_col]])
  if (any(na_weights)) {
    # Fallback: equal weight per row within scenario for unmatched models
    for (sc in unique(df[[scenario_col]][na_weights])) {
      idx_sc <- which(df[[scenario_col]] == sc)
      idx_na <- idx_sc[is.na(df[[weight_col]][idx_sc])]
      if (length(idx_na) > 0) {
        df[[weight_col]][idx_na] <- (1 / length(idx_sc))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Renormalize within each scenario (numerical guard)
  # ---------------------------------------------------------------------------
  df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenario_col))) |>
    dplyr::mutate(
      !!weight_col := .data[[weight_col]] / sum(.data[[weight_col]], na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  # Restore original row order
  df <- df[order(df$.row_id), , drop = FALSE]
  df$.row_id <- NULL

  # ---------------------------------------------------------------------------
  # Validation and diagnostics
  # ---------------------------------------------------------------------------
  if (verbose) {
    weight_sums <- tapply(df[[weight_col]], df[[scenario_col]], sum, na.rm = TRUE)
    bad_sums <- names(weight_sums)[abs(weight_sums - 1) > 1e-6]

    if (length(bad_sums) > 0) {
      warning(
        "Weights do not sum to 1 in scenarios: ",
        paste(bad_sums, collapse = ", "),
        call. = FALSE
      )
    }

    scenarios <- unique(df[[scenario_col]])
    for (sc in scenarios) {
      df_sc <- df[df[[scenario_col]] == sc, , drop = FALSE]
      n_inst <- length(unique(df_sc[[institution_col]][!is.na(df_sc[[institution_col]])]))
      n_mod <- length(unique(df_sc[[model_col]]))
      message(sprintf("[INSTITUTION] %s: %d institutions | %d models", sc, n_inst, n_mod))
    }
  }

  df
}


# =============================================================================
# compute_gcm_weights_by_genealogy
# =============================================================================

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
#' @param clean_col Character. Output column name for cleaned model names. Default "model_clean".
#' @param family_col Character. Output column name for family assignment. Default "model_family".
#' @param weight_col Character. Output column name for weights. Default "w_genealogy".
#' @param keep_original_model Logical. If TRUE, keep original model string in "model_raw".
#' @param verbose Logical. If TRUE, prints mapping/coverage summary and validates weights.
#'
#' @return
#' `gcm_data` with added columns:
#' - `clean_col`: cleaned CMIP model name used for mapping
#' - `family_col`: genealogy family (from Kuma table; fallback to model_clean if unmatched)
#' - `weight_col`: per-row weight within scenario_col; sums to 1 per scenario
#'
#' @examples
#' \dontrun{
#' result <- compute_gcm_weights_by_genealogy(
#'   gcm_data = my_gcm_data,
#'   kuma_table = read.csv("kuma_table_s1.csv"),
#'   method = "family"
#' )
#' }
#'
#' @export
compute_gcm_weights_by_genealogy <- function(
    gcm_data,
    kuma_table,
    model_col = "model",
    scenario_col = "scenario",
    cmip_phase = c("CMIP6", "CMIP5", "CMIP3"),
    method = c("family", "family_sqrt", "independence"),
    clean_col = "model_clean",
    family_col = "model_family",
    weight_col = "w_genealogy",
    keep_original_model = TRUE,
    verbose = TRUE
) {

  cmip_phase <- match.arg(cmip_phase)
  method <- match.arg(method)

  # ---------------------------------------------------------------------------
  # Internal helpers
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
    if (is.na(x) || !nzchar(x)) return(character(0))
    out <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
    out <- trimws(out)
    out[nzchar(out)]
  }

  .build_mapping <- function(kuma_table, cmip_phase, family_col_out) {
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

    map <- data.frame(
      cmip_name = cmip_name,
      stringsAsFactors = FALSE
    )
    map[[family_col_out]] <- family

    # Deduplicate (keep first occurrence)
    map <- map[!duplicated(map$cmip_name), , drop = FALSE]

    # Defensive check: ensure no duplicates remain
    if (any(duplicated(map$cmip_name))) {
      stop("Internal error: duplicate cmip_name entries in mapping table.", call. = FALSE)
    }

    map
  }

  # ---------------------------------------------------------------------------
  # Vectorized weight computation helpers
  # ---------------------------------------------------------------------------

  .compute_weights_vectorized <- function(df, scenario_col, clean_col, family_col,
                                          weight_col, method) {
    # Step 1: Build model -> family mapping per scenario (unique models)
    model_fam_map <- df |>
      dplyr::select(dplyr::all_of(c(scenario_col, clean_col, family_col))) |>
      dplyr::distinct()

    # Step 2: Count models per family per scenario
    fam_model_counts <- model_fam_map |>
      dplyr::group_by(
        dplyr::across(dplyr::all_of(c(scenario_col, family_col)))
      ) |>
      dplyr::summarise(
        .n_models_in_fam = dplyr::n(),
        .groups = "drop"
      )

    # Step 3: Count families per scenario
    fam_counts <- fam_model_counts |>
      dplyr::group_by(dplyr::across(dplyr::all_of(scenario_col))) |>
      dplyr::summarise(
        .n_fam = dplyr::n(),
        .groups = "drop"
      )

    # Step 4: Compute family weights based on method
    if (method == "family") {
      # Equal weight per family
      fam_weights <- fam_model_counts |>
        dplyr::left_join(fam_counts, by = scenario_col) |>
        dplyr::mutate(
          .w_family = 1 / .data$.n_fam
        )
    } else if (method == "family_sqrt") {
      # Weight proportional to 1/sqrt(n_models_in_family)
      fam_weights <- fam_model_counts |>
        dplyr::left_join(fam_counts, by = scenario_col) |>
        dplyr::mutate(
          .w_family_raw = 1 / sqrt(.data$.n_models_in_fam)
        ) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(scenario_col))) |>
        dplyr::mutate(
          .w_family = .data$.w_family_raw / sum(.data$.w_family_raw)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-".w_family_raw")
    } else {
      # independence: equal weight per model (family weight = n_models_in_fam / total_models)
      total_models <- model_fam_map |>
        dplyr::group_by(dplyr::across(dplyr::all_of(scenario_col))) |>
        dplyr::summarise(.n_total = dplyr::n(), .groups = "drop")

      fam_weights <- fam_model_counts |>
        dplyr::left_join(total_models, by = scenario_col) |>
        dplyr::mutate(
          .w_family = .data$.n_models_in_fam / .data$.n_total
        ) |>
        dplyr::select(-".n_total")
    }

    # Step 5: Compute model weights (family weight / n_models_in_family)
    model_weights <- model_fam_map |>
      dplyr::left_join(
        fam_weights |> dplyr::select(dplyr::all_of(c(scenario_col, family_col)),
                                     ".n_models_in_fam", ".w_family"),
        by = c(scenario_col, family_col)
      ) |>
      dplyr::mutate(
        .w_model = .data$.w_family / .data$.n_models_in_fam
      ) |>
      dplyr::select(dplyr::all_of(c(scenario_col, clean_col, ".w_model")))

    # Step 6: Join back to original data and split across rows
    df <- df |>
      dplyr::left_join(model_weights, by = c(scenario_col, clean_col)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(scenario_col, clean_col)))) |>
      dplyr::mutate(
        !!weight_col := .data$.w_model / dplyr::n()
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-".w_model")

    # Step 7: Renormalize within scenario (numerical guard)
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(scenario_col))) |>
      dplyr::mutate(
        !!weight_col := .data[[weight_col]] / sum(.data[[weight_col]], na.rm = TRUE)
      ) |>
      dplyr::ungroup()

    df
  }

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (!is.data.frame(gcm_data)) {
    stop("'gcm_data' must be a data.frame.", call. = FALSE)
  }

  .stop_if_missing(gcm_data, c(model_col, scenario_col), "gcm_data")

  if (nrow(gcm_data) == 0) {
    warning("Input gcm_data has zero rows.", call. = FALSE)
    gcm_data[[clean_col]] <- character(0)
    gcm_data[[family_col]] <- character(0)
    gcm_data[[weight_col]] <- numeric(0)
    return(gcm_data)
  }

  # ---------------------------------------------------------------------------
  # Build mapping and attach families
  # ---------------------------------------------------------------------------

  map <- .build_mapping(kuma_table, cmip_phase = cmip_phase, family_col_out = family_col)

  # Work on a copy with row ID for stable ordering
  df <- gcm_data
  df$.row_id <- seq_len(nrow(df))

  if (keep_original_model && model_col != "model_raw") {
    df$model_raw <- df[[model_col]]
  }

  # Clean model names
  df[[clean_col]] <- .clean_model_name(df[[model_col]])

  # Attach family via left join (safer than merge for preserving order)
  df <- df |>
    dplyr::left_join(
      map,
      by = stats::setNames("cmip_name", clean_col)
    )

  # Fallback for unmatched models: treat as own family
  unmatched <- is.na(df[[family_col]]) | !nzchar(df[[family_col]])
  if (any(unmatched)) {
    df[[family_col]][unmatched] <- df[[clean_col]][unmatched]
  }

  # ---------------------------------------------------------------------------
  # Compute weights (vectorized)
  # ---------------------------------------------------------------------------

  df <- .compute_weights_vectorized(
    df = df,
    scenario_col = scenario_col,
    clean_col = clean_col,
    family_col = family_col,
    weight_col = weight_col,
    method = method
  )

  # Restore original row order
  df <- df[order(df$.row_id), , drop = FALSE]
  df$.row_id <- NULL

  # ---------------------------------------------------------------------------
  # Diagnostics and validation
  # ---------------------------------------------------------------------------

  if (verbose) {
    # Mapping summary
    n_total_models <- length(unique(df[[clean_col]]))
    n_mapped <- sum(unique(df[[clean_col]]) %in% map$cmip_name)
    n_unmapped <- n_total_models - n_mapped

    message(sprintf(
      "[GENEALOGY] CMIP phase=%s | method=%s | unique models=%d | mapped=%d | unmapped=%d (fallback family=model)",
      cmip_phase, method, n_total_models, n_mapped, n_unmapped
    ))

    # Per-scenario diagnostics
    scenarios <- sort(unique(as.character(df[[scenario_col]])))
    for (sc in scenarios) {
      df_sc <- df[df[[scenario_col]] == sc, , drop = FALSE]
      n_fam_sc <- length(unique(df_sc[[family_col]]))
      n_mod_sc <- length(unique(df_sc[[clean_col]]))
      message(sprintf("  %s: %d families | %d models", sc, n_fam_sc, n_mod_sc))
    }

    # Weight sum validation
    weight_sums <- tapply(df[[weight_col]], df[[scenario_col]], sum, na.rm = TRUE)
    bad_sums <- names(weight_sums)[abs(weight_sums - 1) > 1e-6]

    if (length(bad_sums) > 0) {
      warning(
        "Weights do not sum to 1 in scenarios: ",
        paste(bad_sums, collapse = ", "),
        call. = FALSE
      )
    }
  }

  df
}

