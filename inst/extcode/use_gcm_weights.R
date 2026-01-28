# =============================================================================
# Test Script: Verify Improved GCM Weighting Functions
# =============================================================================

library(dplyr)

cat("=======================================================================\n")
cat("Testing Improved GCM Weighting Functions\n")
cat("=======================================================================\n\n")

# =============================================================================
# 1. Create Test Data
# =============================================================================

set.seed(123)

# Simulated GCM ensemble with realistic structure
gcm_data <- data.frame(
  model = c(
    # Family A (Australian): 2 models, ACCESS-CM2 has 2 runs
    "ACCESS-CM2", "ACCESS-CM2", "ACCESS-ESM1-5",
    # Family B (GFDL): 2 models
    "GFDL-CM4", "GFDL-ESM4",
    # Family C (IPSL): 1 model, 3 runs
    "IPSL-CM6A-LR", "IPSL-CM6A-LR", "IPSL-CM6A-LR",
    # Family D (MPI): 2 models
    "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
    # Family E (UK): 2 models
    "UKESM1-0-LL", "HadGEM3-GC31-LL"
  ),
  institution = c(
    "CSIRO", "CSIRO", "CSIRO",
    "NOAA-GFDL", "NOAA-GFDL",
    "IPSL", "IPSL", "IPSL",
    "MPI-M", "MPI-M",
    "MOHC", "MOHC"
  ),
  variant = c(
    "r1i1p1f1", "r2i1p1f1", "r1i1p1f1",
    "r1i1p1f1", "r1i1p1f1",
    "r1i1p1f1", "r2i1p1f1", "r3i1p1f1",
    "r1i1p1f1", "r1i1p1f1",
    "r1i1p1f1", "r1i1p1f1"
  ),
  scenario = rep("SSP2-4.5", 12),
  delta_T = rnorm(12, 2.5, 0.5),
  stringsAsFactors = FALSE
)

# Add a second scenario (subset of models)
gcm_data_ssp5 <- gcm_data[c(1, 3, 5, 6, 9, 11), ]
gcm_data_ssp5$scenario <- "SSP5-8.5"
gcm_data_ssp5$delta_T <- rnorm(6, 4.5, 0.7)

gcm_data <- rbind(gcm_data, gcm_data_ssp5)

cat("=== Test Data Summary ===\n")
cat("Total rows:", nrow(gcm_data), "\n")
cat("Scenarios:", paste(unique(gcm_data$scenario), collapse = ", "), "\n")
cat("Unique models:", length(unique(gcm_data$model)), "\n\n")

# =============================================================================
# 2. Create Mock Kuma Table
# =============================================================================

kuma_table <- data.frame(
  Model = c("ACCESS", "GFDL", "IPSL", "MPI", "UKMO"),
  Family = c("A", "B", "C", "D", "E"),
  `CMIP6 names` = c(
    "ACCESS-CM2, ACCESS-ESM1-5",
    "GFDL-CM4, GFDL-ESM4",
    "IPSL-CM6A-LR",
    "MPI-ESM1-2-HR, MPI-ESM1-2-LR",
    "UKESM1-0-LL, HadGEM3-GC31-LL"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# =============================================================================
# Test 1: Institution-Based Weighting
# =============================================================================

cat("=== Test 1: Institution-Based Weighting ===\n\n")

result_inst <- compute_gcm_weights_by_institution(
  data = gcm_data,
  scenario_col = "scenario",
  institution_col = "institution",
  model_col = "model",
  weight_col = "w_inst",
  verbose = TRUE
)

cat("\nWeight sums:", tapply(result_inst$w_inst, result_inst$scenario, sum), "\n")

cat("\nModel weights (SSP2-4.5):\n")
result_inst %>%
  filter(scenario == "SSP2-4.5") %>%
  group_by(institution, model) %>%
  summarise(w = sum(w_inst), .groups = "drop") %>%
  arrange(institution) %>%
  print(n = 20)

# =============================================================================
# Test 2: Genealogy-Based Weighting (Family)
# =============================================================================

cat("\n=== Test 2: Genealogy-Based Weighting (family) ===\n\n")

result_gen <- compute_gcm_weights_by_genealogy(
  gcm_data = gcm_data,
  kuma_table = kuma_table,
  method = "family",
  verbose = TRUE
)

cat("\nWeight sums:", tapply(result_gen$w_genealogy, result_gen$scenario, sum), "\n")

cat("\nFamily weights (SSP2-4.5):\n")
result_gen %>%
  filter(scenario == "SSP2-4.5") %>%
  group_by(model_family) %>%
  summarise(w = sum(w_genealogy), n_models = n_distinct(model_clean), .groups = "drop") %>%
  print()

# =============================================================================
# Test 3: Genealogy-Based Weighting (Family-Sqrt)
# =============================================================================

cat("\n=== Test 3: Genealogy-Based Weighting (family_sqrt) ===\n\n")

result_sqrt <- compute_gcm_weights_by_genealogy(
  gcm_data = gcm_data,
  kuma_table = kuma_table,
  method = "family_sqrt",
  verbose = TRUE
)

cat("\nFamily weights (SSP2-4.5) - sqrt method:\n")
result_sqrt %>%
  filter(scenario == "SSP2-4.5") %>%
  group_by(model_family) %>%
  summarise(w = sum(w_genealogy), n_models = n_distinct(model_clean), .groups = "drop") %>%
  print()

# =============================================================================
# Test 4: Custom Column Names
# =============================================================================

cat("\n=== Test 4: Custom Column Names ===\n\n")

result_custom <- compute_gcm_weights_by_genealogy(
  gcm_data = gcm_data,
  kuma_table = kuma_table,
  clean_col = "gcm_name",
  family_col = "gcm_family",
  weight_col = "weight",
  verbose = FALSE
)

cat("Columns:", paste(names(result_custom), collapse = ", "), "\n")
cat("Has 'gcm_name':", "gcm_name" %in% names(result_custom), "\n")
cat("Has 'gcm_family':", "gcm_family" %in% names(result_custom), "\n")
cat("Has 'weight':", "weight" %in% names(result_custom), "\n")

# =============================================================================
# Test 5: Weight Validation Utility
# =============================================================================

cat("\n=== Test 5: Weight Validation ===\n\n")

validate_weight_sums(result_gen, "scenario", "w_genealogy")

# =============================================================================
# Test 6: Performance (Large Ensemble)
# =============================================================================

cat("\n=== Test 6: Performance (Large Ensemble) ===\n\n")

n_models <- 50
n_variants <- 5
n_scenarios <- 4

large_data <- expand.grid(
  model_id = 1:n_models,
  variant_id = 1:n_variants,
  scenario = paste0("SSP", 1:n_scenarios)
)
large_data$model <- paste0("MODEL-", large_data$model_id)
large_data$institution <- paste0("INST-", (large_data$model_id %% 10) + 1)

large_kuma <- data.frame(
  Model = paste0("M", 1:n_models),
  Family = paste0("FAM-", (1:n_models %% 15) + 1),
  `CMIP6 names` = paste0("MODEL-", 1:n_models),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Large dataset:", nrow(large_data), "rows\n\n")

t1 <- system.time({
  r1 <- compute_gcm_weights_by_institution(large_data, verbose = FALSE)
})
cat("Institution weighting:", round(t1["elapsed"], 3), "sec\n")

t2 <- system.time({
  r2 <- compute_gcm_weights_by_genealogy(large_data, large_kuma, verbose = FALSE)
})
cat("Genealogy weighting:", round(t2["elapsed"], 3), "sec\n")

# Verify
cat("\nValidation:\n")
validate_weight_sums(r2, "scenario", "w_genealogy")

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("All tests completed!\n")
cat("=======================================================================\n")
