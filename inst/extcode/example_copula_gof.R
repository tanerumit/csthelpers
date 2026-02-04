
# ==============================================================================
# Test: Verify GoF Implementation
# ==============================================================================

test_copula_gof <- function() {

  set.seed(42)

  cat("=== Test 1: Well-specified Gaussian copula ===\n")

  # Generate data from true Gaussian copula
  n <- 200
  rho_true <- 0.6

  # Bivariate normal -> uniform marginals
  z1 <- rnorm(n)
  z2 <- rho_true * z1 + sqrt(1 - rho_true^2) * rnorm(n)
  u_obs <- pnorm(z1)
  v_obs <- pnorm(z2)

  gof1 <- copula_goodness_of_fit(u_obs, v_obs, rho = rho_true)

  cat(sprintf("  RMSE: %.4f (should be small, ~0.01-0.03)\n", gof1$rmse))
  cat(sprintf("  Poor fit: %s\n", gof1$poor_fit))
  cat("\n")

  cat("=== Test 2: Misspecified rho ===\n")

  gof2 <- copula_goodness_of_fit(u_obs, v_obs, rho = 0.0)  # Wrong rho

  cat(sprintf("  RMSE: %.4f (should be larger)\n", gof2$rmse))
  cat(sprintf("  Poor fit: %s\n", gof2$poor_fit))
  if (!is.null(gof2$poor_fit_reason)) cat(sprintf("  Reason: %s\n", gof2$poor_fit_reason))
  cat("\n")


  cat("=== Test 3: Data with tail dependence (Clayton-like) ===\n")

  # Simulate Clayton copula (lower tail dependence) via conditional method
  # Clayton copula: C(u,v) = (u^-theta + v^-theta - 1)^(-1/theta)
  theta <- 3  # Strong lower tail dependence

  u_clayton <- runif(n)
  # Conditional v|u for Clayton
  v_clayton <- (u_clayton^(-theta) * (runif(n)^(-theta/(theta+1)) - 1) + 1)^(-1/theta)

  # Fit Gaussian copula (wrong model)
  z1_c <- qnorm(pmax(1e-10, pmin(1-1e-10, u_clayton)))
  z2_c <- qnorm(pmax(1e-10, pmin(1-1e-10, v_clayton)))
  rho_fit <- cor(z1_c, z2_c)

  gof3 <- copula_goodness_of_fit(u_clayton, v_clayton, rho = rho_fit)

  cat(sprintf("  Fitted rho: %.3f\n", rho_fit))
  cat(sprintf("  Overall RMSE: %.4f\n", gof3$rmse))
  cat(sprintf("  Lower tail RMSE: %.4f\n", gof3$tail_lower_rmse))
  cat(sprintf("  Upper tail RMSE: %.4f\n", gof3$tail_upper_rmse))
  cat(sprintf("  Poor fit: %s\n", gof3$poor_fit))
  if (!is.null(gof3$poor_fit_reason)) cat(sprintf("  Reason: %s\n", gof3$poor_fit_reason))
  cat("\n")


  cat("=== Test 4: Independence ===\n")

  u_indep <- runif(n)
  v_indep <- runif(n)

  gof4 <- copula_goodness_of_fit(u_indep, v_indep, rho = 0)

  cat(sprintf("  RMSE: %.4f (should be small)\n", gof4$rmse))
  cat(sprintf("  Poor fit: %s\n", gof4$poor_fit))
  cat("\n")


  cat("=== Test 5: BVN CDF accuracy check ===\n")

  # Compare against known values
  # P(Z1 < 0, Z2 < 0; rho=0.5) should be ~0.333
  p_test <- .scenwgt_bvn_cdf_scalar(0, 0, 0.5)
  cat(sprintf("  P(Z1<0, Z2<0 | rho=0.5) = %.4f (expected ~0.333)\n", p_test))

  # P(Z1 < 1, Z2 < 1; rho=0) = pnorm(1)^2 = 0.841^2 = 0.708
  p_test2 <- .scenwgt_bvn_cdf_scalar(1, 1, 0)
  expected2 <- pnorm(1)^2
  cat(sprintf("  P(Z1<1, Z2<1 | rho=0) = %.4f (expected %.4f)\n", p_test2, expected2))

  cat("\n=== All tests completed ===\n")

  invisible(list(gof1 = gof1, gof2 = gof2, gof3 = gof3, gof4 = gof4))
}

# Run tests
results <- test_copula_gof()
