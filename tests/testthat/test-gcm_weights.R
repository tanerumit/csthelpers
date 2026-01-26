# Functions: compute_gcm_weights_by_institution (R/gcm_weights.R), compute_gcm_weights_by_genealogy (R/gcm_weights.R)

testthat::test_that("compute_gcm_weights_by_institution splits by institution and model rows", {
  df <- data.frame(
    scenario = c("S1", "S1", "S1", "S2"),
    institution = c("A", "A", "B", NA),
    model = c("m1", "m1", "m2", "m3"),
    stringsAsFactors = FALSE
  )

  out <- compute_gcm_weights_by_institution(df)

  s1 <- out[out$scenario == "S1", "w_inst"]
  testthat::expect_equal(sum(s1), 1)
  testthat::expect_equal(s1, c(0.25, 0.25, 0.5), tolerance = 1e-12)

  s2 <- out[out$scenario == "S2", "w_inst"]
  testthat::expect_equal(s2, 1)
})

testthat::test_that("compute_gcm_weights_by_genealogy returns weights per scenario and keeps model_raw", {
  gcm_data <- data.frame(
    scenario = rep(c("S1", "S2"), each = 3),
    model = c("X1", "X2", "X1", "X1", "X3", "X4"),
    stringsAsFactors = FALSE
  )

  kuma_table <- data.frame(
    Model = c("M1", "M2"),
    Family = c("FamA", "FamB"),
    `CMIP6 names` = c("X1, X2", "X3"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out <- compute_gcm_weights_by_genealogy(
    gcm_data = gcm_data,
    kuma_table = kuma_table,
    cmip_phase = "CMIP6",
    method = "family",
    keep_original_model = TRUE,
    verbose = FALSE
  )

  testthat::expect_true(all(c("model_clean", "model_family", "w_genealogy", "model_raw") %in% names(out)))
  testthat::expect_equal(sum(out$w_genealogy[out$scenario == "S1"]), 1, tolerance = 1e-12)
  testthat::expect_equal(sum(out$w_genealogy[out$scenario == "S2"]), 1, tolerance = 1e-12)
})
