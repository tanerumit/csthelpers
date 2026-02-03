# Functions: estimate_scenario_probs_kde (R/scenario_probabilities.R)

testthat::test_that("area_weight regular yields resolution-invariant masses", {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("dplyr")

  set.seed(42)
  ensemble_data <- data.frame(
    scenario = rep("S1", 60),
    tavg = rnorm(60, mean = 1.5, sd = 0.4),
    prcp = rnorm(60, mean = -0.5, sd = 0.7)
  )

  grid_coarse <- expand.grid(
    tavg = seq(-0.5, 3.5, by = 1),
    prcp = seq(-3, 2, by = 1)
  )

  grid_fine <- expand.grid(
    tavg = seq(-0.5, 3.5, by = 0.5),
    prcp = seq(-3, 2, by = 0.5)
  )

  w_coarse <- estimate_scenario_probs_kde(
    ensemble_data = ensemble_data,
    scenario_grid = grid_coarse,
    group_col = "scenario",
    ta_col = "tavg",
    pr_col = "prcp",
    bw = c(0.6, 0.9),
    normalize = TRUE,
    area_weight = "regular",
    verbose = FALSE
  )

  w_fine <- estimate_scenario_probs_kde(
    ensemble_data = ensemble_data,
    scenario_grid = grid_fine,
    group_col = "scenario",
    ta_col = "tavg",
    pr_col = "prcp",
    bw = c(0.6, 0.9),
    normalize = TRUE,
    area_weight = "regular",
    verbose = FALSE
  )

  nearest_value <- function(x, ref) {
    ref[vapply(x, function(v) which.min(abs(ref - v)), 1L)]
  }

  coarse_ta <- sort(unique(grid_coarse$tavg))
  coarse_pr <- sort(unique(grid_coarse$prcp))

  fine_map <- data.frame(
    tavg = nearest_value(grid_fine$tavg, coarse_ta),
    prcp = nearest_value(grid_fine$prcp, coarse_pr),
    weight = w_fine$S1
  )

  fine_agg <- stats::aggregate(weight ~ tavg + prcp, data = fine_map, sum)
  coarse_df <- data.frame(
    tavg = grid_coarse$tavg,
    prcp = grid_coarse$prcp,
    weight = w_coarse$S1
  )

  merged <- merge(coarse_df, fine_agg, by = c("tavg", "prcp"), sort = TRUE)
  testthat::expect_equal(merged$weight.x, merged$weight.y, tolerance = 5e-2)
})
