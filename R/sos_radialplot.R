#' Create a facetted radial (spider) plot with fixed tick placement
#'
#' @description
#' `sos_radialplot()` generates a facetted radial (polar) bar plot for visualizing
#' categorical groups arranged around a circle (e.g., hydrological indicators, sites, or flow
#' components) and their corresponding numeric values under different scenarios or panels.
#' 
#' The function ensures that:
#' * The radial tick marks (e.g., 10, 20, 30, 40) are consistently placed at the same
#'   angular position (“GAP” sector) across all panels.
#' * The “GAP” sector is hidden from the plot (no spoke or bar drawn).
#' * Tick labels are fixed at a vertical (90°) position by default.
#'
#' @param data A data frame containing the variables to plot.
#' @param axis_col Variable name mapped to the circular axis (e.g., `"group"` or `"site"`).
#' @param panel_col Variable name used to facet panels (e.g., `"scenario"`).
#' @param panel_var Optional variable for nested paneling (currently unused).
#' @param value_col Numeric variable specifying the radial magnitude or bar length.
#' @param cat_colors A named character vector defining the color mapping for `color_col` levels.
#'   Example: `c("No" = "green", "Low" = "yellow", "Medium" = "orange", "High" = "red")`.
#' @param color_col Variable defining color categories for bars.
#' @param radial_breaks Optional numeric vector specifying fixed breaks for the radial axis
#'   (e.g., `seq(0, 50, 10)`). If `NULL`, breaks are automatically computed.
#' @param start_angle Numeric value (in degrees) defining the angular position of the “GAP”
#'   sector where radial tick labels are drawn. Default is `90` (vertical at top).
#'
#' @details
#' The function internally adds a dummy `"GAP"` level to the circular axis to reserve a fixed
#' angular position for radial tick labels. This guarantees consistent tick placement and
#' orientation across multiple facets, regardless of the number of variables.
#'
#' The result is especially useful for comparing hydrological indicators or risk classes
#' across different climate or management scenarios.
#'
#' @return
#' A `ggplot` object (using polar coordinates) that can be further customized with
#' additional `ggplot2` layers or themes.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#'
#' # Example data
#' df <- data.frame(
#'   group = rep(c("HF", "MF", "LF", "RFC", "IF", "Overall"), times = 4),
#'   scenario = rep(c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"), each = 6),
#'   value = c(5, 10, 25, 0, 5, 10,
#'             8, 12, 30, 10, 0, 15,
#'             10, 18, 35, 15, 5, 20,
#'             15, 20, 40, 18, 10, 25),
#'   risk_class = rep(c("Low", "Low", "Medium", "No", "Low", "Low"), times = 4)
#' )
#'
#' cat_colors <- c("No"="green", "Low"="yellow", "Medium"="orange", "High"="red")
#'
#' sos_radialplot(
#'   data = df,
#'   axis_col = group,
#'   panel_col = scenario,
#'   value_col = value,
#'   color_col = risk_class,
#'   cat_colors = cat_colors,
#'   radial_breaks = seq(0, 50, 10),
#'   start_angle = 90
#' )
#' }
#'
#' @import ggplot2 dplyr tibble
#' @export
sos_radialplot <- function(
    data,
    axis_col,           # around-the-circle variable (e.g., group or site)
    panel_col,          # facetting variable (e.g., scenario)
    panel_var = NULL,   # optional (not used)
    value_col,          # radial variable (numeric)
    cat_colors,         # named vector of fill colors
    color_col,          # variable defining color category
    radial_breaks = NULL,   # optional: set fixed radial tick values (e.g., seq(0,50,10))
    start_angle = 90        # angle where GAP (tick labels) will be placed, degrees
) {
  
  # Required packages
  require(ggplot2)
  require(dplyr)
  require(tibble)
  
  # Plotting variables
  lab_size <- 3
  
  # --- Prepare data ---
  df <- data %>% rename(axis_col  := {{axis_col}}, panel_col := {{panel_col}}, 
         value_col := {{value_col}}, color_col := {{color_col}})
  
  df[["axis_col"]]  <- factor(df[["axis_col"]],
         levels = unique(df[["axis_col"]]), ordered = TRUE)
  df[["panel_col"]] <- factor(df[["panel_col"]],
         levels = unique(df[["panel_col"]]), ordered = TRUE)
  df[["color_col"]] <- factor(df[["color_col"]],
         levels = names(cat_colors), ordered = TRUE)
  
  axis_levels  <- levels(df[["axis_col"]])
  panel_levels <- levels(df[["panel_col"]])
  
  # --- Add GAP row per panel (for tick placement only) ---
  gap_row <- df %>%
    group_by(panel_col) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(axis_col = "GAP", value_col = 0)
  
  df <- bind_rows(df, gap_row)
  df[["axis_col"]] <- factor(df[["axis_col"]],
       levels = c(axis_levels, "GAP"), ordered = TRUE)
  
  # --- Label data for axes around rim ---
  lab_df <- tibble(axis_col = factor(axis_levels, levels = c(axis_levels, "GAP")))
  # --- Radial tick labels at the GAP sector ---
  tick_df <- expand_grid(panel_col = panel_levels,
         axis_col  = factor("GAP", levels = c(axis_levels, "GAP")),
         y = radial_breaks[-1], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) 

  # --- Theme ---
  ggtheme <- theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.y = element_blank(),
      axis.title  = element_text(size = 12, color = "black"),
      panel.spacing = grid::unit(1, "lines"),
      plot.title.position = "plot"
    )
  
  # --- Radial scale ---
  if (is.null(radial_breaks)) {
    maxv <- max(df$value_col, na.rm = TRUE)
    radial_breaks <- pretty(c(0, maxv))
  }
  ymax_ext <- max(radial_breaks)
  
  # --- Build plot ---
  p <- ggplot(df) +
    ggtheme +
    # Dashed helper lines
    geom_segment(data = df %>% filter(axis_col != "GAP"),
      aes(x = axis_col, xend = axis_col, y = 0, yend = ymax_ext),
      linetype = "dashed", color = "gray12") +
    # Bars
    geom_col(data = df %>% filter(axis_col != "GAP"),
      aes(x = axis_col, y = value_col, fill = color_col),
      width = 0.8, color = "gray20", alpha = .9, show.legend = TRUE) +
    # scales
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_continuous(limits = c(0, ymax_ext), breaks = radial_breaks) +
    scale_fill_manual(values = cat_colors, name = "Levels", drop = FALSE) +
    # Other
    facet_wrap(~ panel_col, nrow = 2) +
    coord_polar(clip = "on", start = start_angle * pi / 180) +
    labs(x = NULL, y = NULL) +
    # Labeling
    geom_label(data = lab_df,
        aes(x = axis_col, y = ymax_ext, label = axis_col),
        fill = "white", border.color = "white", color = "black", alpha = 0.8,
        fontface = "bold", size = lab_size, vjust = 0.5) +
    geom_text(data = tick_df, aes(x = axis_col, y = y, label = y),
        size = lab_size, color = "black", alpha = 1, inherit.aes = FALSE)
  
  return(p)
}
