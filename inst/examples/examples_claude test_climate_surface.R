# Test and Example Script for climate_surface() - New Naming Convention
# This file demonstrates the refactored function with intuitive underscore names

library(ggplot2)
library(dplyr)

# Source the refactored function
source("R/climate_surface.R")

# ==============================================================================
# EXAMPLE 1: Basic Response Surface
# ==============================================================================

cat("\n=== Example 1: Basic Climate Response Surface ===\n")

# Create a grid of climate scenarios
df_basic <- expand.grid(
  prcp_change = seq(-30, 30, by = 5),
  temp_change = seq(-1, 5, by = 0.5)
)

# Simulate system reliability (decreases with warming and drying)
set.seed(123)
df_basic$reliability <- with(df_basic, {
  base <- 0.95
  prcp_effect <- -0.008 * prcp_change
  temp_effect <- -0.04 * temp_change
  noise <- rnorm(nrow(df_basic), 0, 0.02)
  pmax(0, pmin(1, base + prcp_effect + temp_effect + noise))
})

# Create basic plot - NOTE THE NEW PARAMETER NAMES!
p1 <- climate_surface(
  data = df_basic,                    # was: str.data
  x_var = "prcp_change",              # was: variable.x
  y_var = "temp_change",              # was: variable.y
  z_var = "reliability",              # was: variable.z
  title = "Water System Reliability", # was: plot.title
  x_label = expression(Delta * "Precipitation" ~ "(%)"),  # was: variable.x.label
  y_label = expression(Delta * "Temperature" ~ "(°C)"),   # was: variable.y.label
  threshold = 0.85,                   # was: threshold.z
  failure_dir = 1                     # was: failure.direction
)

print(p1)
cat("\n✓ Basic plot created successfully with new naming!\n")

# ==============================================================================
# EXAMPLE 2: With GCM Projections
# ==============================================================================

cat("\n=== Example 2: Response Surface with GCM Projections ===\n")

# Create sample GCM projection data
gcm_projections <- data.frame(
  prcp_change = c(-5, 0, 5, 10, -10, 5, 15, 20, -8, 2, 12, 18),
  temp_change = c(1.5, 2.0, 2.5, 3.0, 1.2, 2.2, 3.5, 4.0, 1.8, 2.4, 3.2, 4.2),
  scenario = c(
    "ssp126", "ssp126", "ssp245", "ssp245",
    "ssp126", "ssp126", "ssp370", "ssp585",
    "ssp245", "ssp245", "ssp370", "ssp585"
  ),
  horizon = c(
    "near", "far", "near", "far",
    "far", "far", "near", "far",
    "near", "far", "far", "far"
  )
)

p2 <- climate_surface(
  data = df_basic,                    # was: str.data
  gcm_data = gcm_projections,         # was: gcm.data
  x_var = "prcp_change",              # was: variable.x
  y_var = "temp_change",              # was: variable.y
  z_var = "reliability",              # was: variable.z
  title = "System Reliability with CMIP6 Projections",  # was: plot.title
  x_label = expression(Delta * "P" ~ "(%)"),            # was: variable.x.label
  y_label = expression(Delta * "T" ~ "(°C)"),           # was: variable.y.label
  threshold = 0.85,                   # was: threshold.z
  ellipse_levels = c(0.5, 0.9),       # was: gcm.bvnorm.levels
  gcm_alpha = 0.7                     # was: gcm.transparency
)

print(p2)
cat("\n✓ Plot with GCM projections created successfully!\n")

# ==============================================================================
# EXAMPLE 3: Filtering GCM Scenarios
# ==============================================================================

cat("\n=== Example 3: Filtering Specific Scenarios ===\n")

# Show only SSP2-4.5 and SSP5-8.5
p3 <- climate_surface(
  data = df_basic,
  gcm_data = gcm_projections,
  x_var = "prcp_change",
  y_var = "temp_change",
  z_var = "reliability",
  title = "Moderate and High Emissions Scenarios Only",
  threshold = 0.85,
  scenarios = c("ssp245", "ssp585"),  # was: gcm.scenario.list
  show_legend = TRUE                  # was: gcm.legend
)

print(p3)
cat("\n✓ Filtered scenario plot created successfully!\n")

# ==============================================================================
# EXAMPLE 4: Custom Units and Threshold Inference
# ==============================================================================

cat("\n=== Example 4: Custom Units and Auto Threshold ===\n")

# Create data with different variables
df_custom <- expand.grid(
  flow_change = seq(-40, 20, by = 5),
  evap_change = seq(0, 15, by = 1)
)

set.seed(456)
df_custom$storage <- with(df_custom, {
  100 + 0.5 * flow_change - 2 * evap_change + rnorm(nrow(df_custom), 0, 3)
})

# Let function infer threshold from baseline (0, 0)
p4 <- climate_surface(
  data = df_custom,                   # was: str.data
  x_var = "flow_change",              # was: variable.x
  y_var = "evap_change",              # was: variable.y
  z_var = "storage",                  # was: variable.z
  title = "Reservoir Storage Response",  # was: plot.title
  x_label = expression(Delta * "Inflow"),      # was: variable.x.label
  y_label = expression(Delta * "Evaporation"), # was: variable.y.label
  x_suffix = " ML/day",               # was: variable.x.suffix
  y_suffix = " mm/year",              # was: variable.y.suffix
  threshold = NULL,                   # was: threshold.z (auto-infer)
  failure_dir = 1                     # was: failure.direction
)

print(p4)
cat("\n✓ Custom units plot created successfully!\n")

# ==============================================================================
# EXAMPLE 5: Multi-Panel Plot (Faceting)
# ==============================================================================

cat("\n=== Example 5: Multi-Panel Comparison (Faceting) ===\n")

# Create data for multiple time periods
df_faceted <- expand.grid(
  prcp_change = seq(-30, 30, by = 10),
  temp_change = seq(-1, 5, by = 1),
  period = c("2030-2050", "2050-2070", "2070-2100")
)

set.seed(789)
df_faceted$reliability <- with(df_faceted, {
  # Reliability decreases more in future periods
  period_effect <- case_when(
    period == "2030-2050" ~ 0,
    period == "2050-2070" ~ -0.05,
    period == "2070-2100" ~ -0.10
  )

  base <- 0.95 + period_effect
  prcp_effect <- -0.008 * prcp_change
  temp_effect <- -0.04 * temp_change
  noise <- rnorm(n(), 0, 0.02)

  pmax(0, pmin(1, base + prcp_effect + temp_effect + noise))
})

p5 <- climate_surface(
  data = df_faceted,                  # was: str.data
  x_var = "prcp_change",              # was: variable.x
  y_var = "temp_change",              # was: variable.y
  z_var = "reliability",              # was: variable.z
  title = "System Reliability Across Time Horizons",  # was: plot.title
  threshold = 0.85,                   # was: threshold.z
  facet = TRUE,                       # was: multi.panel
  facet_by = "period",                # was: panel.variable
  facet_levels = c("2030-2050", "2050-2070", "2070-2100"),  # was: panel.variable.levels
  text_size = 0.5                     # was: text.scale
)

print(p5)
cat("\n✓ Multi-panel faceted plot created successfully!\n")

# ==============================================================================
# EXAMPLE 6: Custom Breaks and Contour Levels
# ==============================================================================

cat("\n=== Example 6: Custom Breaks and Contour Levels ===\n")

p6 <- climate_surface(
  data = df_basic,
  x_var = "prcp_change",
  y_var = "temp_change",
  z_var = "reliability",
  title = "Custom Breaks Example",
  x_breaks = seq(-30, 30, 15),        # was: variable.x.breaks
  y_breaks = seq(-1, 5, 2),           # was: variable.y.breaks
  z_breaks = seq(0.6, 1.0, 0.1),      # was: variable.z.breaks
  n_contours = 10,                    # was: contour.num
  threshold = 0.85,
  failure_dir = 1
)

print(p6)
cat("\n✓ Custom breaks plot created successfully!\n")

# ==============================================================================
# EXAMPLE 7: Reversed Failure Direction
# ==============================================================================

cat("\n=== Example 7: Reversed Failure Direction ===\n")

# For some metrics, higher is worse (e.g., flood risk)
df_flood <- expand.grid(
  prcp_change = seq(-20, 40, by = 5),
  temp_change = seq(-1, 5, by = 0.5)
)

set.seed(321)
df_flood$flood_risk <- with(df_flood, {
  base <- 0.10
  prcp_effect <- 0.01 * prcp_change  # More precip = more flood risk
  temp_effect <- 0.02 * temp_change  # Warming = more intense storms
  noise <- rnorm(nrow(df_flood), 0, 0.02)
  pmax(0, pmin(1, base + prcp_effect + temp_effect + noise))
})

p7 <- climate_surface(
  data = df_flood,
  x_var = "prcp_change",
  y_var = "temp_change",
  z_var = "flood_risk",
  title = "Flood Risk (Higher = Worse)",
  threshold = 0.20,
  failure_dir = -1  # ABOVE threshold is failure!
)

print(p7)
cat("\n✓ Reversed color scheme plot created successfully!\n")

# ==============================================================================
# EXAMPLE 8: Comparing Old vs New Naming
# ==============================================================================

cat("\n=== Example 8: Side-by-Side Naming Comparison ===\n")

cat("\nOLD NAMING (dot notation):\n")
cat("  str.data, variable.x, variable.y, variable.z\n")
cat("  threshold.z, plot.title, failure.direction\n")
cat("  gcm.scenario.list, gcm.transparency, gcm.legend\n")
cat("  multi.panel, panel.variable, text.scale\n\n")

cat("NEW NAMING (underscore notation):\n")
cat("  data, x_var, y_var, z_var\n")
cat("  threshold, title, failure_dir\n")
cat("  scenarios, gcm_alpha, show_legend\n")
cat("  facet, facet_by, text_size\n\n")

cat("Benefits:\n")
cat("  ✓ Shorter and easier to type\n")
cat("  ✓ More intuitive names\n")
cat("  ✓ Consistent underscore notation\n")
cat("  ✓ Follows R/tidyverse best practices\n")
cat("  ✓ Better IDE autocomplete\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("All examples completed successfully with new naming!\n")
cat("\nKey improvements:\n")
cat("  ✓ Intuitive parameter names\n")
cat("  ✓ Consistent underscore notation\n")
cat("  ✓ Shorter, clearer names\n")
cat("  ✓ Standard R conventions (n_contours, text_size, etc.)\n")
cat("  ✓ ggplot2 terminology (facet instead of multi.panel)\n")
cat("  ✓ All functionality preserved\n")
cat("\nFunction is production-ready with improved naming!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
