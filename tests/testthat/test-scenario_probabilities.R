# Functions: estimate_scenario_probs_kde (R/scenario_probabilities.R), estimate_scenario_probs_mvnorm (R/scenario_probabilities.R)

testthat::test_that("estimate_scenario_probs_kde normalizes and respects support mask", {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("dplyr")

  set.seed(1)
  ensemble_data <- data.frame(
    scenario = rep(c("S1", "S2"), each = 20),
    tavg = c(rnorm(20, 0, 0.3), rnorm(20, 1, 0.3)),
    prcp = c(rnorm(20, 0, 0.4), rnorm(20, 1, 0.4))
  )

  scenario_grid <- expand.grid(
    tavg = seq(-1, 2, by = 0.5),
    prcp = seq(-1, 2, by = 0.5)
  )

  out <- estimate_scenario_probs_kde(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    bw = c(0.5, 0.5),
    normalize = TRUE,
    area_weight = "regular",
    support = list(ta = c(-0.5, 1.5), pr = c(-0.5, 1.5)),
    verbose = FALSE
  )

  testthat::expect_true(all(c("S1", "S2") %in% names(out)))
  testthat::expect_equal(sum(out$S1), 1, tolerance = 1e-12)
  testthat::expect_equal(sum(out$S2), 1, tolerance = 1e-12)
  outside <- scenario_grid$tavg < -0.5 | scenario_grid$tavg > 1.5 |
    scenario_grid$prcp < -0.5 | scenario_grid$prcp > 1.5
  testthat::expect_true(all(out$S1[outside] == 0))
  testthat::expect_true(all(out$S2[outside] == 0))
})

testthat::test_that("estimate_scenario_probs_mvnorm returns diagnostics when requested", {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("dplyr")

  set.seed(2)
  ensemble_data <- data.frame(
    scenario = rep("S1", 30),
    tavg = rnorm(30, 0.5, 0.2),
    prcp = rnorm(30, -0.2, 0.3)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-0.5, 1.5, by = 0.25),
    prcp = seq(-1, 1, by = 0.25)
  )

  out <- estimate_scenario_probs_mvnorm(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    normalize = TRUE,
    diagnostics = TRUE,
    verbose = FALSE
  )

  testthat::expect_true("S1" %in% names(out))
  testthat::expect_equal(sum(out$S1), 1, tolerance = 1e-12)
  testthat::expect_true(!is.null(attr(out, "mvn_params")))
  testthat::expect_true(!is.null(attr(out, "scaling_params")))
  testthat::expect_true(!is.null(attr(out, "fit_method")))
})
