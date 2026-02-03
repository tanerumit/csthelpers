test_that("copula: basic output structure and normalization", {
  set.seed(21)

  ensemble_data <- data.frame(
    tavg = c(rnorm(90, 1.0, 0.9), rnorm(90, 2.3, 0.8)),
    prcp = c(rnorm(90, 10.0, 2.0), rnorm(90, 13.0, 1.7)),
    scenario = rep(c("SSP1", "SSP2"), each = 90)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-1, 5, by = 0.25),
    prcp = seq(5, 19, by = 0.5)
  )

  out <- compute_scenario_surface_weights_copula(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    ta_col = "tavg",
    pr_col = "prcp",
    group_col = "scenario",
    bw_method = "nrd0",
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c("tavg", "prcp") %in% names(out)))
  expect_true(all(c("SSP1", "SSP2") %in% names(out)))

  expect_equal(sum(out$SSP1), 1, tolerance = 1e-10)
  expect_equal(sum(out$SSP2), 1, tolerance = 1e-10)
})

test_that("copula: diagnostics attributes exist and are consistent", {
  set.seed(22)

  ensemble_data <- data.frame(
    tavg = c(rnorm(70, 0.7, 0.9), rnorm(70, 2.4, 0.8)),
    prcp = c(rnorm(70, 9.2, 1.9), rnorm(70, 13.2, 1.6)),
    scenario = rep(c("SSP1", "SSP2"), each = 70)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-2, 5, by = 0.5),
    prcp = seq(6, 17, by = 0.5)
  )

  out <- compute_scenario_surface_weights_copula(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    bw_method = "nrd0",
    normalize = TRUE,
    area_weight = "none",
    scale = "global",
    diagnostics = TRUE,
    verbose = FALSE
  )

  expect_true(is.character(attr(out, "skipped_groups")))
  expect_true(is.list(attr(out, "bandwidth_used")))
  expect_true(is.list(attr(out, "bandwidth_method")))
  expect_true(is.list(attr(out, "effective_sample_size")))
  expect_true(is.list(attr(out, "scaling_params")))
  expect_true(is.list(attr(out, "copula_rho")))

  rhos <- attr(out, "copula_rho")
  expect_true(all(c("SSP1", "SSP2") %in% names(rhos)))
  expect_true(all(is.finite(unlist(rhos))))
  expect_true(all(abs(unlist(rhos)) < 1))
})

test_that("copula: support masking zeros weights outside bounds", {
  set.seed(23)

  ensemble_data <- data.frame(
    tavg = rnorm(140, 1.4, 1.0),
    prcp = rnorm(140, 11.2, 2.0),
    scenario = rep("SSP1", 140)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-4, 6, by = 0.5),
    prcp = seq(2, 22, by = 1.0)
  )

  support <- list(ta = c(0, 3), pr = c(8, 14))

  out <- compute_scenario_surface_weights_copula(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    support = support,
    normalize = FALSE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  inside <- scenario_grid$tavg >= 0 & scenario_grid$tavg <= 3 &
            scenario_grid$prcp >= 8 & scenario_grid$prcp <= 14

  expect_true(all(out$SSP1[!inside] == 0))
  expect_true(any(out$SSP1[inside] > 0))
})

test_that("copula: weights_col changes output relative to unweighted", {
  set.seed(24)

  n <- 160
  ensemble_data <- data.frame(
    tavg = rnorm(n, 1.1, 1.0),
    prcp = rnorm(n, 10.5, 2.0),
    scenario = rep("SSP1", n)
  )

  # Upweight a warm/wet cluster
  w <- rep(0.3, n)
  idx <- seq_len(30)
  ensemble_data$tavg[idx] <- ensemble_data$tavg[idx] + 2.0
  ensemble_data$prcp[idx] <- ensemble_data$prcp[idx] + 2.5
  w[idx] <- 8
  ensemble_data$w <- w

  scenario_grid <- expand.grid(
    tavg = seq(-2, 6, by = 0.25),
    prcp = seq(5, 18, by = 0.5)
  )

  out_unw <- compute_scenario_surface_weights_copula(
    ensemble_data = subset(ensemble_data, select = -w),
    scenario_grid = scenario_grid,
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  out_w <- compute_scenario_surface_weights_copula(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    weights_col = "w",
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = FALSE,
    verbose = FALSE
  )

  expect_gt(sum(abs(out_unw$SSP1 - out_w$SSP1)), 0.1)
  expect_equal(sum(out_unw$SSP1), 1, tolerance = 1e-10)
  expect_equal(sum(out_w$SSP1), 1, tolerance = 1e-10)
})

test_that("copula: groups are skipped when insufficient samples or near-zero variance", {
  set.seed(25)

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
    compute_scenario_surface_weights_copula(
      ensemble_data = ensemble_data,
      scenario_grid = scenario_grid,
      min_samples = 5L,
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

  expect_true("OK" %in% names(out))
  expect_false("ZEROVAR" %in% names(out))
  expect_false("TOOSMALL" %in% names(out))
})

test_that("copula: input validation errors for missing columns and empty surface", {
  ensemble_data <- data.frame(tavg = 1:5, prcp = 1:5, scenario = "SSP1")
  scenario_grid <- data.frame(tavg = numeric(0), prcp = numeric(0))

  expect_error(
    compute_scenario_surface_weights_copula(
      ensemble_data = ensemble_data,
      scenario_grid = scenario_grid,
      verbose = FALSE
    ),
    "cannot be empty"
  )

  expect_error(
    compute_scenario_surface_weights_copula(
      ensemble_data = data.frame(tavg = 1:5, scenario = "SSP1"),
      scenario_grid = data.frame(tavg = 1, prcp = 1),
      verbose = FALSE
    ),
    "Missing required columns"
  )
})

test_that("copula: rho is near zero for independent marginals (approx)", {
  set.seed(26)

  n <- 200
  ensemble_data <- data.frame(
    tavg = rnorm(n, 0, 1),
    prcp = rnorm(n, 10, 2),  # independent of tavg by construction
    scenario = rep("SSP1", n)
  )

  scenario_grid <- expand.grid(
    tavg = seq(-3, 3, by = 0.5),
    prcp = seq(6, 14, by = 0.5)
  )

  out <- compute_scenario_surface_weights_copula(
    ensemble_data = ensemble_data,
    scenario_grid = scenario_grid,
    bw_method = "nrd0",
    normalize = TRUE,
    area_weight = "none",
    scale = "none",
    diagnostics = TRUE,
    verbose = FALSE
  )

  rho <- attr(out, "copula_rho")$SSP1
  expect_true(is.finite(rho))
  expect_lt(abs(rho), 0.3)  # loose tolerance; KDE-CDF pseudo obs is approximate
})
