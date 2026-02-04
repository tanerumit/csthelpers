#' Compute Weighted Impacts Across Climate Scenarios
#'
#' @description
#' Computes weighted impact summaries by joining scenario weights to impact data
#' and aggregating across realizations.
#'
#' @param scenario_weights Data frame combining scenario realizations with their
#'   climate-based weights. Must contain columns for scenario ID, realization,
#'   climate variables (ta, pr), scenario grouping, and weight.
#' @param impact_data Data frame in long format with impact values. Must contain
#'   scenario ID and value columns, optionally a location column.
#' @param scenario_cols Named list specifying columns in `scenario_weights`:
#'   \describe{
#'     \item{scenario_id}{Links to impact_data (default: "strid")}
#'     \item{realization}{Realization/ensemble member (default: "rlz")}
#'     \item{scenario}{Scenario grouping column(s), NULL for single scenario,
#'       or character vector for multiple (default: "scenario")}
#'     \item{weight}{Weight column (default: "weight")}
#'   }
#' @param impact_cols Named list specifying columns in `impact_data`:
#'   \describe{
#'     \item{scenario_id}{Links to scenario_weights (default: "strid")}
#'     \item{location}{Location/site column, or NULL for single location (default: NULL)}
#'     \item{value}{Impact value column (default: "value")}
#'   }
#' @param realization_agg_fn Function to aggregate across realizations (default: median).
#' @param return_detail Logical; if TRUE, also returns per-realization results.
#'
#' @return Data frame with weighted impacts per scenario (and location if provided).
#'   If `return_detail = TRUE`, returns a list with `summary` and `by_realization`.
#'
#' @examples
#' \dontrun{
#' # Basic usage with defaults
#' result <- compute_weighted_impacts(
#'   scenario_weights = my_scenario_weights,
#'   impact_data = my_impacts
#' )
#'
#' # Custom column names
#' result <- compute_weighted_impacts(
#'   scenario_weights = my_scenario_weights,
#'   impact_data = my_impacts,
#'   scenario_cols = list(
#'     scenario_id = "run_id",
#'     realization = "member",
#'     scenario = c("ssp", "gcm"),
#'     weight = "w"
#'   ),
#'   impact_cols = list(
#'     scenario_id = "run_id",
#'     location = "gauge",
#'     value = "discharge"
#'   )
#' )
#' }
#'
#' @export
compute_weighted_impacts <- function(
    scenario_weights,
    impact_data,
    scenario_cols = list(
      scenario_id = "strid",
      realization = "rlz",
      scenario = "scenario",
      weight = "weight"
    ),
    impact_cols = list(
      scenario_id = "strid",
      location = NULL,
      value = "value"
    ),
    realization_agg_fn = stats::median,
    return_detail = FALSE
) {



  # ---------------------------------------------------------------------------

  # Resolve column specifications with defaults
  # ---------------------------------------------------------------------------

  .fill_defaults <- function(user_list, defaults) {
    for (nm in names(defaults)) {
      if (!nm %in% names(user_list)) {
        user_list[[nm]] <- defaults[[nm]]
      }
    }
    user_list
  }

  scenario_defaults <- list(scenario_id = "strid", realization = "rlz",
                            scenario = "scenario", weight = "weight")
  impact_defaults <- list(scenario_id = "strid", location = NULL, value = "value")

  scenario_cols <- .fill_defaults(scenario_cols, scenario_defaults)
  impact_cols <- .fill_defaults(impact_cols, impact_defaults)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (!is.data.frame(scenario_weights)) {
    stop("scenario_weights must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(impact_data)) {
    stop("impact_data must be a data.frame.", call. = FALSE)
  }

  # Validate scenario_weights columns
  has_scenario <- !is.null(scenario_cols$scenario)
  req_scn <- c(scenario_cols$scenario_id, scenario_cols$realization, scenario_cols$weight)
  if (has_scenario) req_scn <- c(req_scn, scenario_cols$scenario)
  miss_scn <- setdiff(req_scn, names(scenario_weights))
  if (length(miss_scn) > 0L) {
    stop("scenario_weights missing columns (check scenario_cols): ",
         paste(miss_scn, collapse = ", "), call. = FALSE)
  }

  # Validate impact_data columns
  has_location <- !is.null(impact_cols$location)
  req_imp <- c(impact_cols$scenario_id, impact_cols$value)
  if (has_location) req_imp <- c(req_imp, impact_cols$location)
  miss_imp <- setdiff(req_imp, names(impact_data))
  if (length(miss_imp) > 0L) {
    stop("impact_data missing columns (check impact_cols): ",
         paste(miss_imp, collapse = ", "), call. = FALSE)
  }

  if (!is.function(realization_agg_fn)) {
    stop("realization_agg_fn must be a function.", call. = FALSE)
  }

  if (!is.numeric(scenario_weights[[scenario_cols$weight]])) {
    stop("scenario_weights weight column must be numeric.", call. = FALSE)
  }

  if (!is.numeric(impact_data[[impact_cols$value]])) {
    stop("impact_data value column must be numeric.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Prepare scenario weights
  # ---------------------------------------------------------------------------

  if (has_scenario) {
    scn_col <- scenario_cols$scenario
    if (length(scn_col) == 1L) {
      scn_df <- scenario_weights |>
        dplyr::transmute(
          .scenario_id = .data[[scenario_cols$scenario_id]],
          .realization = .data[[scenario_cols$realization]],
          .scenario_key = as.character(.data[[scn_col]]),
          .weight = .data[[scenario_cols$weight]]
        )
    } else {
      # Multiple scenario columns: composite key
      scn_df <- scenario_weights |>
        dplyr::mutate(
          .scenario_key = do.call(paste, c(lapply(scn_col, function(col) .data[[col]]), sep = "|"))
        ) |>
        dplyr::transmute(
          .scenario_id = .data[[scenario_cols$scenario_id]],
          .realization = .data[[scenario_cols$realization]],
          .scenario_key = .data$.scenario_key,
          .weight = .data[[scenario_cols$weight]]
        )
      # Lookup for restoring columns later
      scenario_lookup <- scenario_weights |>
        dplyr::mutate(
          .scenario_key = do.call(paste, c(lapply(scn_col, function(col) .data[[col]]), sep = "|"))
        ) |>
        dplyr::select(dplyr::all_of(c(".scenario_key", scn_col))) |>
        dplyr::distinct()
    }
  } else {
    scn_df <- scenario_weights |>
      dplyr::transmute(
        .scenario_id = .data[[scenario_cols$scenario_id]],
        .realization = .data[[scenario_cols$realization]],
        .scenario_key = ".all",
        .weight = .data[[scenario_cols$weight]]
      )
  }

  # ---------------------------------------------------------------------------
  # Prepare impact data
  # ---------------------------------------------------------------------------

  if (has_location) {
    impact_df <- impact_data |>
      dplyr::transmute(
        .scenario_id = .data[[impact_cols$scenario_id]],
        .location = as.character(.data[[impact_cols$location]]),
        .value = .data[[impact_cols$value]]
      )
  } else {
    impact_df <- impact_data |>
      dplyr::transmute(
        .scenario_id = .data[[impact_cols$scenario_id]],
        .location = ".all",
        .value = .data[[impact_cols$value]]
      )
  }

  # ---------------------------------------------------------------------------
  # Join scenario_weights to impact_data
  # ---------------------------------------------------------------------------

  joined <- dplyr::inner_join(scn_df, impact_df, by = ".scenario_id")

  if (nrow(joined) == 0L) {
    stop(
      "No rows after joining. Check that scenario_id values match ",
      "between scenario_weights and impact_data.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Compute weighted impacts
  # ---------------------------------------------------------------------------

  by_realization <- joined |>
    dplyr::mutate(.weighted_value = .data$.weight * .data$.value) |>
    dplyr::group_by(.data$.scenario_key, .data$.realization, .data$.location) |>
    dplyr::summarise(
      weighted_impact = sum(.data$.weighted_value, na.rm = TRUE),
      .groups = "drop"
    )

  summary_df <- by_realization |>
    dplyr::group_by(.data$.scenario_key, .data$.location) |>
    dplyr::summarise(
      weighted_impact = realization_agg_fn(.data$weighted_impact, na.rm = TRUE),
      n_realizations = dplyr::n(),
      .groups = "drop"
    )

  # ---------------------------------------------------------------------------
  # Restore original column names
  # ---------------------------------------------------------------------------

  # Scenario columns
  if (!has_scenario) {
    summary_df <- summary_df |> dplyr::select(-".scenario_key")
    by_realization <- by_realization |> dplyr::select(-".scenario_key")
  } else if (length(scenario_cols$scenario) == 1L) {
    summary_df <- summary_df |>
      dplyr::rename(!!scenario_cols$scenario := .data$.scenario_key)
    by_realization <- by_realization |>
      dplyr::rename(!!scenario_cols$scenario := .data$.scenario_key)
  } else {
    summary_df <- summary_df |>
      dplyr::left_join(scenario_lookup, by = ".scenario_key") |>
      dplyr::select(-".scenario_key")
    by_realization <- by_realization |>
      dplyr::left_join(scenario_lookup, by = ".scenario_key") |>
      dplyr::select(-".scenario_key")
  }

  # Location column
  if (!has_location) {
    summary_df <- summary_df |> dplyr::select(-".location")
    by_realization <- by_realization |>
      dplyr::rename(!!scenario_cols$realization := .data$.realization) |>
      dplyr::select(-".location")
  } else {
    summary_df <- summary_df |>
      dplyr::rename(!!impact_cols$location := .data$.location)
    by_realization <- by_realization |>
      dplyr::rename(
        !!scenario_cols$realization := .data$.realization,
        !!impact_cols$location := .data$.location
      )
  }

  summary_df <- summary_df |> dplyr::arrange(dplyr::across(dplyr::everything()))

  if (isTRUE(return_detail)) {
    return(list(summary = summary_df, by_realization = by_realization))
  }

  summary_df
}
