#' Calculate Environmental Flow Risk Assessment (ERFA) Class
#'
#' @description
#' Compares the 70 eco-hydrological indicators (from `calc_eflow_metrics()`)
#' between two discharge scenarios (e.g., baseline and future) and classifies
#' the environmental flow (E-flow) risk for each indicator group and overall.
#'
#' @param baseline_tbl A data frame returned by [calc_eflow_metrics()] representing
#' the baseline (reference) discharge scenario.
#' @param scenario_tbl A data frame returned by [calc_eflow_metrics()] representing
#' the scenario to compare against the baseline.
#' @param sensi_threshold Numeric; sensitivity threshold for percent change (default = 30).
#' Indicators with an absolute percent change greater than this threshold
#' are considered “different” (sensitive).
#'
#' @details
#' The function calculates percent change for each indicator:
#' \deqn{pct\_change = 100 * (value_{scenario} - value_{baseline}) / value_{baseline}}
#'
#' Indicators whose absolute percent change exceeds the sensitivity threshold
#' are counted as “different”. The counts are summarized by indicator group
#' (HF, MF, LF, RFC, IF) and overall. Based on the number of indicators outside
#' the threshold range, the function assigns an environmental flow risk class.
#'
#' ### Risk classification thresholds
#'
#' | group | 0 (No change) | 1 (Low) | 2 (Medium) | 3 (High) |
#' |--------|---------------|---------|-------------|-----------|
#' | **Overall** | 0 | 1–22 | 23–46 | 47–70 |
#' | **HF** | 0 | 1–5 | 6–10 | 11–16 |
#' | **MF** | 0 | 1–8 | 9–16 | 17–24 |
#' | **LF** | 0 | 1–5 | 6–10 | 11–16 |
#' | **RFC** | 0 | 1–2 | 3–5 | 6–8 |
#' | **IF** | 0 | 1–2 | 3–4 | 5–6 |
#'
#' Risk class codes and suggested colors:
#' \itemize{
#'   \item 0 = No change — 🟦 Blue
#'   \item 1 = Low risk — 🟩 Green
#'   \item 2 = Medium risk — 🟧 Amber
#'   \item 3 = High risk — 🟥 Red
#' }
#'
#' @return
#' A list with two data frames:
#' \describe{
#'   \item{detail}{
#'     Per-indicator comparison table with baseline and scenario values,
#'     percent change, and logical flag (`outside_range`) indicating if
#'     the change exceeds the threshold.
#'   }
#'   \item{summary}{
#'     Summary table by group and overall, containing:
#'     \itemize{
#'       \item `group` — indicator group (HF, MF, LF, RFC, IF, Overall)
#'       \item `n_total` — total indicators in the group
#'       \item `n_outside` — count of indicators outside threshold
#'       \item `perc_outside` — percent of indicators outside threshold
#'       \item `risk_class` — numeric class code (0–3)
#'       \item `risk_label` — descriptive label with color code
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' dates <- seq.Date(as.Date("2000-01-01"), as.Date("2019-12-31"), by = "day")
#' Q_base <- pmax(0, 40 + 15*sin(2*pi*(1:length(dates))/365) + rnorm(length(dates),0,5))
#' Q_fut  <- pmax(0, 42 + 17*sin(2*pi*(1:length(dates))/365) + rnorm(length(dates),0,6))
#'
#' base_tbl <- calc_eflow_metrics(Q_base, dates)
#' scen_tbl <- calc_eflow_metrics(Q_fut,  dates)
#'
#' erfa <- calculate_erfa_class(base_tbl, scen_tbl, sensi_threshold = 30)
#' erfa$summary
#' }
#'
#' @import dplyr
#' @export
calculate_erfa_class <- function(baseline_tbl,
                                 scenario_tbl,
                                 sensi_threshold = 30) {
  
  # --- Merge and compute percent change
  df <- baseline_tbl %>%
    select(name, group, value_base = value) %>%
    inner_join(scenario_tbl %>% select(name, value_scen = value),
               by = "name") %>%
    mutate(
      pct_change = 100 * (value_scen - value_base) / value_base,
      outside_range = abs(pct_change) > sensi_threshold
    )
  
  # --- Count indicators outside sensitivity range by group
  summary_tbl <- df %>%
    group_by(group) %>%
    summarise(
      n_total = n(),
      n_outside = sum(outside_range, na.rm = TRUE),
      perc_outside = 100 * n_outside / n_total,
      .groups = "drop"
    )
  
  # --- Add overall summary
  total_row <- summary_tbl %>%
    summarise(
      group = "Overall",
      n_total = sum(n_total),
      n_outside = sum(n_outside),
      perc_outside = 100 * n_outside / n_total
    )
  
  summary_tbl <- bind_rows(summary_tbl, total_row)
  
  # --- Define risk classification thresholds
  risk_thresholds <- tribble(
    ~group, ~low, ~medium, ~high,
    "Overall", 1, 23, 47,
    "HF", 1, 6, 11,
    "MF", 1, 9, 17,
    "LF", 1, 6, 11,
    "RFC", 1, 3, 6,
    "IF", 1, 3, 5
  )
  
  # --- Assign risk class and labels
  summary_tbl <- summary_tbl %>%
    left_join(risk_thresholds, by = "group") %>%
    mutate(
      risk_class = case_when(
        n_outside == 0 ~ 0,
        n_outside >= high ~ 3,
        n_outside >= medium ~ 2,
        n_outside >= low ~ 1,
        TRUE ~ NA_real_
      ),
      risk_label = factor(
        risk_class,
        levels = 0:3,
        labels = c("No change (🟦 Blue)",
                   "Low (🟩 Green)",
                   "Medium (🟧 Amber)",
                   "High (🟥 Red)")
      )
    ) %>%
    select(group, n_total, n_outside, perc_outside, risk_class, risk_label) %>%
    mutate(
      group = factor(group,
                     levels = c("HF", "MF", "LF", "RFC", "IF", "Overall"),
                     ordered = TRUE)
    ) %>%
    arrange(group)
  
  # --- Return both detailed and summary outputs
  return(list(detail = df, summary = summary_tbl))
}
