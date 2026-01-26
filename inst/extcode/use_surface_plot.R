

# Data
str_lookup <- read_csv("data/s4w_rhine_strlookup.csv")
str_out <- read_csv("data/s4w_rhine_strtest_Q50.csv") %>%
  left_join(str_lookup, by = "strid") %>%
  filter(rlz == 2) %>%
  select(strid, tavg, prcp, Lobith)

gcm_data <- read_csv("data/s4w_rhine_gcm_summary_stats.csv")
gcm_use <- gcm_data %>% filter(horizon == "far") %>% na.omit()

data = str_out
x_var = "prcp"
y_var = "tavg"
z_var = "Lobith"
threshold = NULL
title = "Climate Response Surface"
x_label = expression(Delta ~ "Precipitation")
y_label = expression(Delta ~ "Temperature")
x_suffix = "%"
y_suffix = "\u00B0C"
failure_dir = 1
x_breaks = NULL
y_breaks = NULL
n_contours = 25
z_limits = NULL
panel_size_in = 6.0
legend_barheight_in = 0.25
text_size = 0.7
facet = FALSE
facet_by = NULL
facet_levels = NULL

range(str_out$Lobith)



p <- climate_surface_base(
    data = str_out,
    x_var = "prcp",
    y_var = "tavg",
    z_var = "Lobith",
    threshold = 1200,
    title = "Climate Response Surface",
    x_label = expression(Delta ~ "Precipitation"),
    y_label = expression(Delta ~ "Temperature"),
    x_suffix = "%",
    y_suffix = "\u00B0C",
    failure_dir = 1,
    x_breaks = NULL,
    y_breaks = NULL,
    n_contours = 13,
    z_limits = c(600, 3200),
    legend_barwidth_in = p_legend_barwidth,
    legend_barheight_in = p_legend_barheight,
    text_size = text_size,
    facet = FALSE,
    facet_by = NULL,
    facet_levels = NULL)


# Report
p_width <- 6.5
p_height <- 6
p_legend_barwidth <- 4.31
p_legend_barheight = 0.25
text_size = 0.95

# Powerpoint
# p_width <- 6.5 * plot_scaler
# p_height <- 6  * plot_scaler
# p_legend_barwidth <- 4.31  * 2.1
# p_legend_barheight = 0.25  * plot_scaler
# text_size = 0.95 * 1.8

ggsave("TEMP/figure.png", p, width = p_width, height = p_height, dpi = 300)
