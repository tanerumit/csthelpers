# Functions: compute_weighted_impacts (R/scenario_weighted_impacts.R)

testthat::test_that("compute_weighted_impacts aggregates across realizations", {
  testthat::skip_if_not_installed("dplyr")

  scenario_key <- data.frame(
    strid = c(1, 2, 1, 2),
    rlz = c(1, 1, 2, 2),
    tavg = c(0, 1, 0, 1),
    prcp = c(0, 1, 0, 1)
  )

  kde_weights <- data.frame(
    tavg = c(0, 1),
    prcp = c(0, 1),
    scenario = c("S1", "S1"),
    weight = c(0.6, 0.4)
  )

  impacts <- data.frame(
    strid = c(1, 2),
    siteA = c(10, 20)
  )

  out <- suppressWarnings(compute_weighted_impacts(
    scenario_key = scenario_key,
    kde_weights = kde_weights,
    impacts = impacts,
    site_col = "siteA",
    return_detail = TRUE
  ))

  testthat::expect_true(all(c("summary", "by_realization") %in% names(out)))
  testthat::expect_equal(out$summary$weighted_impact, 14, tolerance = 1e-12)
  testthat::expect_equal(out$summary$n_realizations, 2)
})
