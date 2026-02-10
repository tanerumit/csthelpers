# Functions: compute_weighted_impacts (R/weighted_impacts.R)

testthat::test_that("compute_weighted_impacts aggregates across realizations", {
  testthat::skip_if_not_installed("dplyr")

  scenario_weights <- data.frame(
    strid = c(1, 2, 1, 2),
    rlz = c(1, 1, 2, 2),
    scenario = "S1",
    weight = c(0.6, 0.4, 0.6, 0.4)
  )

  impact_data <- data.frame(
    strid = c(1, 2),
    value = c(10, 20)
  )

  out <- suppressWarnings(compute_weighted_impacts(
    scenario_weights = scenario_weights,
    impact_data = impact_data,
    return_detail = TRUE
  ))

  testthat::expect_true(all(c("summary", "by_realization") %in% names(out)))
  testthat::expect_equal(out$summary$weighted_impact, 14, tolerance = 1e-12)
  testthat::expect_equal(out$summary$n_realizations, 2)
})
