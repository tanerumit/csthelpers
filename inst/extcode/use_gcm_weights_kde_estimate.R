

# ------------------------------------------------------------------------------
# Simple end-to-end test: KDE weights with vs without institution weighting
# ------------------------------------------------------------------------------

# Assumes these are already sourced in your session:
# - estimate_scenario_probs_kde()  (the revised function)
# - compute_gcm_weights_by_institution() (from gcm_weights*.R)
#
# If not:

source("R/gcm_weights.R")
source("R/scenario_probabilities.R")
source("R/gcm_weights_utils.R")

set.seed(42)


# GCM-Projections
gcm_data <- read_csv("data/s4w_rhine_gcm_summary_stats.csv") %>% filter(horizon == "far")


n <- 60L
ensemble_data <- data.frame(
  scenario     = rep(c("SSP1", "SSP2"), each = n/2),
  Model        = paste0("GCM_", sprintf("%02d", 1:n)),
  Institution  = c(
    rep("Inst_A", 25),  # over-represented
    rep("Inst_B", 15),
    rep("Inst_C", 10),
    rep("Inst_D", 10)
  ),
  stringsAsFactors = FALSE
)

# Give the two scenarios different clouds (so you can see separation)
ensemble_data$tavg <- ifelse(
  ensemble_data$scenario == "SSP1",
  rnorm(n, mean = 1.6, sd = 0.35),
  rnorm(n, mean = 2.6, sd = 0.45)
)

ensemble_data$prcp <- ifelse(
  ensemble_data$scenario == "SSP1",
  rnorm(n, mean =  4.5, sd = 1.2),
  rnorm(n, mean = -0.5, sd = 1.4)
)

# ------------------------------------------------------------------------------
# 2) Define a stress-test grid
# ------------------------------------------------------------------------------

scenario_grid <- expand.grid(
  tavg = seq(0.5, 3.8, by = 0.25),
  prcp = seq(-5,  8,   by = 0.5)
)

# ------------------------------------------------------------------------------
# 3) Compute institution weights and attach as a per-row weight column
# ------------------------------------------------------------------------------

# You need to pass the appropriate column names expected by your implementation.
# Most commonly: model_col / institution_col (or similar).
# If your function uses different argument names, adjust here.

w_inst <- compute_gcm_weights_by_institution(
  data = ensemble_data,
  model_col          = "Model",
  institution_col    = "Institution"
)

# Expectation: it returns either:
# - a numeric vector aligned with ensemble_data rows, OR
# - a data.frame with a weight column aligned by row / Model.
#
# Handle both patterns:

if (is.numeric(w_inst) && length(w_inst) == nrow(ensemble_data)) {
  ensemble_data$w_inst <- w_inst
} else if (is.data.frame(w_inst)) {

  # Try the most likely weight column name(s)
  w_col <- intersect(names(w_inst), c("weight", "w", "w_inst", "inst_weight", "gcm_weight"))
  if (length(w_col) != 1) {
    stop("Cannot infer weight column name from compute_gcm_weights_by_institution() output.")
  }

  # Try joining by Model if present; else assume row-aligned
  if ("Model" %in% names(w_inst)) {
    ensemble_data <- merge(ensemble_data, w_inst[, c("Model", w_col)], by = "Model", all.x = TRUE, sort = FALSE)
    names(ensemble_data)[names(ensemble_data) == w_col] <- "w_inst"
    # merge() reorders; restore original order
    ensemble_data <- ensemble_data[match(paste0("GCM_", sprintf("%02d", 1:n)), ensemble_data$Model), ]
  } else {
    if (nrow(w_inst) != nrow(ensemble_data)) stop("Weight table not row-aligned and has no Model column.")
    ensemble_data$w_inst <- w_inst[[w_col]]
  }

} else {
  stop("Unsupported output type from compute_gcm_weights_by_institution().")
}

# Sanity checks
stopifnot(is.numeric(ensemble_data$w_inst))
stopifnot(all(is.finite(ensemble_data$w_inst)))
stopifnot(all(ensemble_data$w_inst >= 0))
stopifnot(any(ensemble_data$w_inst > 0))

# ------------------------------------------------------------------------------
# 4) KDE without weights
# ------------------------------------------------------------------------------

w_unweighted <- estimate_scenario_probs_kde(
  ensemble_data = ensemble_data,
  scenario_grid = scenario_grid,
  group_col     = "scenario",
  ta_col        = "tavg",
  pr_col        = "prcp",
  weights_col   = NULL,
  scale         = "global",
  area_weight   = "regular",
  normalize     = TRUE,
  diagnostics   = TRUE,
  verbose       = FALSE
)

# ------------------------------------------------------------------------------
# 5) KDE with institution weights
# ------------------------------------------------------------------------------

w_weighted <- estimate_scenario_probs_kde(
  ensemble_data = ensemble_data,
  scenario_grid = scenario_grid,
  group_col     = "scenario",
  ta_col        = "tavg",
  pr_col        = "prcp",
  weights_col   = "w_inst",
  scale         = "global",
  area_weight   = "regular",
  normalize     = TRUE,
  diagnostics   = TRUE,
  verbose       = FALSE
)

print(tibble(uni = w_unweighted$SSP2, wgt = w_weighted$SSP2), n = 378)

# ------------------------------------------------------------------------------
# 6) Compare results quickly
# ------------------------------------------------------------------------------

cat("\nColumn sums (should be ~1 each when normalize=TRUE):\n")
print(rbind(
  unweighted = colSums(w_unweighted[, c("SSP1", "SSP2")]),
  weighted   = colSums(w_weighted[,   c("SSP1", "SSP2")])
))

cat("\nEffective sample sizes (diagnostics):\n")
print(list(
  unweighted = attr(w_unweighted, "effective_sample_size"),
  weighted   = attr(w_weighted,   "effective_sample_size")
))

cat("\nSurface similarity (correlation across grid points):\n")
cor_ssp1 <- stats::cor(w_unweighted$SSP1, w_weighted$SSP1)
cor_ssp2 <- stats::cor(w_unweighted$SSP2, w_weighted$SSP2)
print(c(SSP1 = cor_ssp1, SSP2 = cor_ssp2))

# Optional: see where the largest changes occur (top 10 grid points)
delta_ssp1 <- abs(w_weighted$SSP1 - w_unweighted$SSP1)
delta_ssp2 <- abs(w_weighted$SSP2 - w_unweighted$SSP2)

top1 <- order(delta_ssp1, decreasing = TRUE)[1:10]
top2 <- order(delta_ssp2, decreasing = TRUE)[1:10]

cat("\nTop 10 grid points by |delta| (SSP1):\n")
print(cbind(scenario_grid[top1, ], delta = delta_ssp1[top1]))

cat("\nTop 10 grid points by |delta| (SSP2):\n")
print(cbind(scenario_grid[top2, ], delta = delta_ssp2[top2]))
