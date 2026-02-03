

# fit surfaces
w_kde <- compute_scenario_surface_weights_kde(...)
w_mvn <- compute_scenario_surface_weights_mvnorm(...)
w_cop <- compute_scenario_surface_weights_copula(...)

# score surfaces (same evaluator)
gof_kde <- evaluate_scenario_surface_weights(w_kde, ensemble_data, mapping="knn", k=5)
gof_mvn <- evaluate_scenario_surface_weights(w_mvn, ensemble_data, mapping="knn", k=5)
gof_cop <- evaluate_scenario_surface_weights(w_cop, ensemble_data, mapping="knn", k=5)

# compare side-by-side externally
rbind(
  transform(gof_kde, method="kde"),
  transform(gof_mvn, method="mvnorm"),
  transform(gof_cop, method="copula")
)
