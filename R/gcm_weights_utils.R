
# =============================================================================
# Utility: Validate weight sums
# =============================================================================

#' Validate that weights sum to 1 within each scenario
#'
#' @param data data.frame with weights
#' @param scenario_col Scenario column name
#' @param weight_col Weight column name
#' @param tol Tolerance for sum check (default 1e-6)
#'
#' @return Logical TRUE if all scenarios sum to 1 (within tolerance), FALSE otherwise.
#'   Also prints diagnostic information.
#'
#' @export
validate_weight_sums <- function(data, scenario_col, weight_col, tol = 1e-6) {
  if (!scenario_col %in% names(data)) {
    stop("Column '", scenario_col, "' not found.", call. = FALSE)
  }
  if (!weight_col %in% names(data)) {
    stop("Column '", weight_col, "' not found.", call. = FALSE)
  }

  weight_sums <- tapply(data[[weight_col]], data[[scenario_col]], sum, na.rm = TRUE)

  cat("Weight sums by scenario:\n")
  for (sc in names(weight_sums)) {
    status <- if (abs(weight_sums[sc] - 1) <= tol) "OK" else "FAIL"
    cat(sprintf("  %s: %.10f [%s]\n", sc, weight_sums[sc], status))
  }

  all_ok <- all(abs(weight_sums - 1) <= tol)

  if (all_ok) {
    cat("All scenarios pass validation.\n")
  } else {
    cat("WARNING: Some scenarios do not sum to 1.\n")
  }

  invisible(all_ok)
}


# =============================================================================
# Utility: Summarize weights by model
# =============================================================================

#' Summarize weights aggregated to model level
#'
#' @param data data.frame with weights
#' @param scenario_col Scenario column name
#' @param model_col Model column name
#' @param weight_col Weight column name
#' @param family_col Optional family column name for grouping
#'
#' @return data.frame with model-level weight summary
#'
#' @export
summarize_model_weights <- function(
    data,
    scenario_col,
    model_col,
    weight_col,
    family_col = NULL
) {
  group_cols <- c(scenario_col, model_col)
  if (!is.null(family_col) && family_col %in% names(data)) {
    group_cols <- c(group_cols, family_col)
  }

  result <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      weight_total = sum(.data[[weight_col]], na.rm = TRUE),
      n_rows = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::across(dplyr::all_of(scenario_col)),
      dplyr::desc(.data$weight_total)
    )

  result
}
