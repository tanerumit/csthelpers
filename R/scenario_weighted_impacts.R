

compute_weighted_impacts <- function(
    scenario_key,
    kde_weights,
    impacts,
    site_col,
    scenario_col = "scenario",
    ta_col = "tavg",
    pr_col = "prcp",
    rlz_col = "rlz",
    strid_col = "strid",
    weight_col = "weight",
    agg_across_rlz = stats::median,
    return_detail = FALSE
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  # ---- dplyr join helper (silence expected many-to-many warning) ----
  .inner_join_safe <- function(x, y, by) {
    ij <- dplyr::inner_join
    fmls <- tryCatch(names(formals(ij)), error = function(e) character(0))

    if ("relationship" %in% fmls) {
      ij(x, y, by = by, relationship = "many-to-many")
    } else {
      suppressWarnings(ij(x, y, by = by))
    }
  }

  # ---- validation ----
  if (!is.data.frame(scenario_key)) stop("scenario_key must be a data.frame.", call. = FALSE)
  if (!is.data.frame(kde_weights)) stop("kde_weights must be a data.frame.", call. = FALSE)
  if (!is.data.frame(impacts)) stop("impacts must be a data.frame.", call. = FALSE)

  if (!is.character(site_col) || length(site_col) != 1 || !nzchar(site_col)) {
    stop("site_col must be a non-empty character scalar.", call. = FALSE)
  }

  req_key <- c(strid_col, rlz_col, ta_col, pr_col)
  miss_key <- setdiff(req_key, names(scenario_key))
  if (length(miss_key) > 0) {
    stop("scenario_key missing required columns: ", paste(miss_key, collapse = ", "), call. = FALSE)
  }

  req_kde <- c(ta_col, pr_col, scenario_col, weight_col)
  miss_kde <- setdiff(req_kde, names(kde_weights))
  if (length(miss_kde) > 0) {
    stop("kde_weights missing required columns: ", paste(miss_kde, collapse = ", "), call. = FALSE)
  }

  req_imp <- c(strid_col, site_col)
  miss_imp <- setdiff(req_imp, names(impacts))
  if (length(miss_imp) > 0) {
    stop("impacts missing required columns: ", paste(miss_imp, collapse = ", "), call. = FALSE)
  }

  if (!is.function(agg_across_rlz)) {
    stop("agg_across_rlz must be a function (e.g., stats::median).", call. = FALSE)
  }

  # ---- ensure KDE weights are unique per (tavg, prcp, scenario) ----
  kde_core <- kde_weights |>
    dplyr::transmute(
      .ta = .data[[ta_col]],
      .pr = .data[[pr_col]],
      .scenario = as.character(.data[[scenario_col]]),
      .w = .data[[weight_col]]
    )

  if (!is.numeric(kde_core$.w)) stop("kde_weights[[weight_col]] must be numeric.", call. = FALSE)

  kde_check <- kde_core |>
    dplyr::group_by(.data$.ta, .data$.pr, .data$.scenario) |>
    dplyr::summarise(
      n_unique_w = dplyr::n_distinct(.data$.w),
      .groups = "drop"
    )

  if (any(kde_check$n_unique_w > 1)) {
    stop(
      "kde_weights has inconsistent weights for the same (tavg, prcp, scenario). ",
      "Provide a single unique set of weights per (tavg, prcp, scenario).",
      call. = FALSE
    )
  }

  kde_distinct <- kde_core |>
    dplyr::group_by(.data$.ta, .data$.pr, .data$.scenario) |>
    dplyr::summarise(weight = dplyr::first(.data$.w), .groups = "drop")

  # ---- join: scenario_key (all rlz) x KDE weights (scenarios) x impacts (by strid) ----
  # This join is intentionally many-to-many on (.ta, .pr): multiple realizations per grid
  # and multiple scenarios per grid point.
  df_key <- scenario_key |>
    dplyr::transmute(
      strid = .data[[strid_col]],
      rlz = .data[[rlz_col]],
      .ta = .data[[ta_col]],
      .pr = .data[[pr_col]]
    )

  df <- .inner_join_safe(
    x = df_key,
    y = kde_distinct,
    by = c(".ta", ".pr")
  ) |>
    dplyr::inner_join(
      impacts |>
        dplyr::select(
          strid = dplyr::all_of(strid_col),
          impact = dplyr::all_of(site_col)
        ),
      by = "strid"
    )

  if (nrow(df) == 0) {
    stop(
      "No rows after joining. Check that (tavg, prcp) in scenario_key matches kde_weights, ",
      "and that strid in scenario_key matches impacts.",
      call. = FALSE
    )
  }

  if (!is.numeric(df$impact)) stop("Selected site impact column must be numeric.", call. = FALSE)

  # ---- compute weighted impact per (scenario, rlz) then aggregate across rlz ----
  by_rlz <- df |>
    dplyr::mutate(weighted = .data$weight * .data$impact) |>
    dplyr::group_by(.data$.scenario, .data$rlz) |>
    dplyr::summarise(
      weighted_impact = sum(.data$weighted, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(!!scenario_col := .data$.scenario)

  summary <- by_rlz |>
    dplyr::group_by(.data[[scenario_col]]) |>
    dplyr::summarise(
      weighted_impact = agg_across_rlz(.data$weighted_impact, na.rm = TRUE),
      n_realizations = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data[[scenario_col]])

  if (isTRUE(return_detail)) {
    return(list(summary = summary, by_realization = by_rlz))
  }

  summary
}
