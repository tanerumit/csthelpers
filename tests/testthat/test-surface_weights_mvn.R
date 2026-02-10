# Functions: compute_scenario_surface_weights_mvn (R/scenario_surface_weights_mvn.R)

test_that("mvnorm: basic output structure and normalization", {
  set.seed(1)

  ensemble_data <- data.frame(
    tavg = c(rnorm(60, 1, 0.8), rnorm(60, 2, 0.7)),
    prcp = c(rnorm(60, 10, 2.0), rnorm(60, 12, 1.8)),
    scenario = rep(c("SSP1", "SSP2"), each = 60)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-1, 5, by = 0.25),
    prcp = seq(5, 18, by = 0.5)
  )

  out <- compute_scenario_surface_weights_mvn(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    ta_col = "tavg",
    pr_col = "prcp",
    group_col = "scenario",
    robust = TRUE,
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c("tavg", "prcp", "scenario", "weight") %in% names(out)))

  sums <- tapply(out$weight, out$scenario, sum)
  expect_true(all(c("SSP1", "SSP2") %in% names(sums)))
  expect_equal(sums[["SSP1"]], 1, tolerance = 1e-10)
  expect_equal(sums[["SSP2"]], 1, tolerance = 1e-10)

  expect_true(all(out$weight >= 0 | is.na(out$weight)))
})

test_that("mvnorm: diagnostics attributes exist and are consistent", {
  set.seed(2)

  ensemble_data <- data.frame(
    tavg = c(rnorm(40, 0.5, 1.0), rnorm(40, 2.2, 0.9)),
    prcp = c(rnorm(40, 9.5, 1.7), rnorm(40, 13.0, 1.6)),
    scenario = rep(c("SSP1", "SSP2"), each = 40)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-2, 5, by = 0.5),
    prcp = seq(6, 17, by = 0.5)
  )

  out <- compute_scenario_surface_weights_mvn(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    robust = TRUE,
    normalize = TRUE,
    area_weight = "none",
    scale = "global",
    diagnostics = TRUE,
    verbose = FALSE
  )

  expect_true(is.character(attr(out, "skipped_groups")))
  expect_true(is.list(attr(out, "mvn_params")))
  expect_true(is.list(attr(out, "scaling_params")))
  expect_true(is.list(attr(out, "fit_method")))
  expect_true(is.list(attr(out, "effective_sample_size")))

  mvn_params <- attr(out, "mvn_params")
  expect_true(all(c("SSP1", "SSP2") %in% names(mvn_params)))

  # mu is length-2, Sigma is 2x2
  expect_equal(length(mvn_params$SSP1$mu), 2)
  expect_equal(dim(mvn_params$SSP1$Sigma), c(2, 2))
})

test_that("mvnorm: support masking zeros weights outside bounds", {
  set.seed(3)

  ensemble_data <- data.frame(
    tavg = rnorm(80, 1.5, 1.0),
    prcp = rnorm(80, 11.0, 2.0),
    scenario = rep("SSP1", 80)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-4, 6, by = 0.5),
    prcp = seq(2, 22, by = 1.0)
  )

  support <- list(ta = c(0, 3), pr = c(8, 14))

  out <- compute_scenario_surface_weights_mvn(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    support = support,
    normalize = FALSE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  inside <- out$tavg >= 0 & out$tavg <= 3 &
    out$prcp >= 8 & out$prcp <= 14 &
    out$scenario == "SSP1"

  expect_true(all(out$weight[!inside & out$scenario == "SSP1"] == 0))
  expect_true(any(out$weight[inside] > 0))
})

test_that("mvnorm: robust=TRUE rejects weights_col", {
  set.seed(4)

  ensemble_data <- data.frame(
    tavg = rnorm(20, 0, 1),
    prcp = rnorm(20, 10, 2),
    scenario = rep("SSP1", 20),
    w = runif(20)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-2, 2, by = 0.5),
    prcp = seq(6, 14, by = 1.0)
  )

  expect_error(
    compute_scenario_surface_weights_mvn(
      ensemble_data = ensemble_data,
      scenario_grid = scenario_grid,
      weights_col = "w",
      robust = TRUE,
      verbose = FALSE
    ),
    "not compatible"
  )
})

test_that("mvnorm: weights_col changes the output when robust=FALSE", {
  set.seed(5)

  n <- 80
  ensemble_data <- data.frame(
    tavg = rnorm(n, 1, 1),
    prcp = rnorm(n, 10, 2),
    scenario = rep("SSP1", n)
  )

  # Strongly upweight a cluster (shifted warm/wet)
  w <- rep(0.2, n)
  idx <- seq_len(20)
  ensemble_data$tavg[idx] <- ensemble_data$tavg[idx] + 2.5
  ensemble_data$prcp[idx] <- ensemble_data$prcp[idx] + 3.0
  w[idx] <- 5
  ensemble_data$w <- w

  scenario_grid <- expand.grid(
    tavg = seq(-2, 6, by = 0.25),
    prcp = seq(5, 18, by = 0.5)
  )

  out_unw <- compute_scenario_surface_weights_mvn(
    ensemble_data = subset(ensemble_data, select = -w),
    scenario_grid = scenario_grid,
    robust = FALSE,
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  out_w <- compute_scenario_surface_weights_mvn(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    weights_col = "w",
    robust = FALSE,
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  # Outputs should differ measurably
  expect_gt(sum(abs(out_unw$weight - out_w$weight)), 0.1)

  # Both should normalize
  expect_equal(sum(out_unw$weight), 1, tolerance = 1e-10)
  expect_equal(sum(out_w$weight), 1, tolerance = 1e-10)
})

test_that("mvnorm: groups are skipped when insufficient samples or near-zero variance", {
  set.seed(6)

  ensemble_data <- data.frame(
    tavg = c(rnorm(10, 0, 1), rep(1, 10), rnorm(3, 2, 0.1)),
    prcp = c(rnorm(10, 10, 2), rnorm(10, 10, 2), rnorm(3, 12, 0.1)),
    scenario = c(rep("OK", 10), rep("ZEROVAR", 10), rep("TOOSMALL", 3))
  )

  scenario_grid <- expand.grid(
    tavg = seq(-3, 5, by = 0.5),
    prcp = seq(6, 16, by = 0.5)
  )

  out <- suppressWarnings(
    compute_scenario_surface_weights_mvn(
      ensemble_data = ensemble_data,
      scenario_grid = scenario_grid,
      min_samples = 5L,
      robust = TRUE,
      normalize = TRUE,
      area_weight = "none",
      scale = "none",
      diagnostics = TRUE,
      verbose = FALSE
    )
  )

  skipped <- attr(out, "skipped_groups")
  expect_true("ZEROVAR" %in% skipped)
  expect_true("TOOSMALL" %in% skipped)
  expect_false("OK" %in% skipped)

  expect_true("OK" %in% out$scenario)
  expect_false("ZEROVAR" %in% out$scenario)
  expect_false("TOOSMALL" %in% out$scenario)
})

test_that("mvnorm: input validation errors for missing columns and empty grid", {
  ensemble_data <- data.frame(tavg = 1:5, prcp = 1:5, scenario = "SSP1")
  scenario_grid <- data.frame(tavg = numeric(0), prcp = numeric(0))

  expect_error(
    compute_scenario_surface_weights_mvn(
      ensemble_data = ensemble_data,
      scenario_grid = scenario_grid,
      verbose = FALSE
    ),
    "cannot be empty"
  )

  expect_error(
    compute_scenario_surface_weights_mvn(
      ensemble_data = data.frame(tavg = 1:5, scenario = "SSP1"),
      scenario_grid = data.frame(tavg = 1, prcp = 1),
      verbose = FALSE
    ),
    "Missing required columns"
  )
})
