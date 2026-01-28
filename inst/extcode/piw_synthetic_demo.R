# =============================================================================
# Synthetic demo: compare GCM weighting approaches
# - Equal (uniform)
# - Institution
# - Genealogy
# - Performance–Independence Weighting (PIW)
#
# This script expects the three scripts to be available locally and sourced:
#   - gcm_weights.R            (genealogy + institution structural weights)
#   - gcm_weights_utils.R      (optional helpers)
#   - gcm_weights_pi.R         (PIW framework; compute_gcm_weights_bma() + stats)
#
# NOTE: gcm_weights.R depends on dplyr and rlang.
# =============================================================================

# --- Paths (edit if needed) ---
path_genealogy_institution <- "R/gcm_weights.R"
path_utils                <- "R/gcm_weights_utils.R"
path_piw                  <- "R/gcm_weights_pi.R"

source(path_genealogy_institution)
source(path_utils)
source(path_piw)

# =============================================================================
# 1) Synthetic data generator
# =============================================================================
make_synthetic_gcm_example <- function(
  seed = 42,
  n_time_hist = 120,          # months
  n_draws_future = 200,       # e.g., gridcells / bootstrap draws per model-scenario
  scenarios_future = c("SSP1-2.6", "SSP2-4.5", "SSP5-8.5")
) {
  set.seed(seed)

  models <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-CM6-1", "CNRM-ESM2-1",
              "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", "GFDL-CM4")

  institutions <- c("CSIRO", "CSIRO", "CNRM", "CNRM", "IPSL", "MIROC", "MPI", "NOAA-GFDL")
  inst_df <- data.frame(model = models, institution = institutions, stringsAsFactors = FALSE)

  # Genealogy table (minimal Kuma-like structure)
  kuma_table <- data.frame(
    Model = paste0("M", seq_along(models)),
    Family = c("FamA","FamA","FamB","FamB","FamC","FamD","FamE","FamF"),
    stringsAsFactors = FALSE
  )

  # IMPORTANT: add CMIP6 names column explicitly, exact spelling
  kuma_table[["CMIP6 names"]] <- models

  # --- Historical obs (truth) ---
  t <- seq_len(n_time_hist)
  obs <- 10 + 2 * sin(2*pi*t/12) + 0.01 * t + rnorm(n_time_hist, sd = 0.6)  # seasonal + trend + noise
  obs_data <- data.frame(time = t, y = obs, stringsAsFactors = FALSE)

  # --- Historical model series (vary bias + noise to create different skill) ---
  # Lower rmse => better. We create a range.
  bias_vec  <- c(0.3, 0.6, -0.2, 0.1, -0.7, 0.0, 0.4, -0.1)
  noise_vec <- c(0.7, 0.8, 0.65, 0.7, 0.9, 0.6, 0.75, 0.7)

  hist_list <- list()
  for (i in seq_along(models)) {
    m <- models[i]
    yhat <- obs + bias_vec[i] + rnorm(n_time_hist, sd = noise_vec[i])
    hist_list[[i]] <- data.frame(
      model = m,
      scenario = "historical",
      time = t,
      y = yhat,
      stringsAsFactors = FALSE
    )
  }
  gcm_hist <- do.call(rbind, hist_list)
  gcm_hist <- merge(gcm_hist, inst_df, by = "model", all.x = TRUE)

  # --- Compute model-level skill (RMSE vs obs) on historical only ---
  # IMPORTANT: compute_skill_metric() expects equal length vectors (we aligned by time index).
  skill_df <- data.frame(model = models, skill = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(models)) {
    m <- models[i]
    mod_vals <- gcm_hist$y[gcm_hist$model == m]
    skill_df$skill[i] <- compute_skill_metric(mod_vals, obs, skill_metric = "rmse")
  }

  # --- Future projections (single variable "proj") ---
  # Create scenario-dependent shifts; some models respond more strongly.
  scenario_shift <- c("SSP1-2.6" = 0.5, "SSP2-4.5" = 1.2, "SSP5-8.5" = 2.4)
  model_sens <- c(0.9, 1.1, 1.0, 1.2, 0.8, 1.3, 0.95, 1.05)

  fut_list <- list()
  k <- 1
  for (sc in scenarios_future) {
    for (i in seq_along(models)) {
      m <- models[i]
      mu <- mean(obs) + scenario_shift[[sc]] * model_sens[i]
      # Draws represent spatial cells / annual means / bootstrap replicates, etc.
      vals <- rnorm(n_draws_future, mean = mu, sd = 0.8)
      fut_list[[k]] <- data.frame(
        model = m,
        scenario = sc,
        draw = seq_len(n_draws_future),
        proj = vals,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  gcm_future <- do.call(rbind, fut_list)
  gcm_future <- merge(gcm_future, inst_df, by = "model", all.x = TRUE)

  list(
    gcm_future = gcm_future,
    gcm_hist = gcm_hist,
    obs_data = obs_data,
    skill_df = skill_df,
    kuma_table = kuma_table
  )
}

# =============================================================================
# 2) Convenience runner: compute 4 weighting options + summary stats
# =============================================================================
compute_all_weighting_options <- function(
  gcm_data,
  skill_data,
  kuma_table,
  model_col = "model",
  scenario_col = "scenario",
  institution_col = "institution",
  skill_col = "skill",
  target_col = "proj",
  # PIW hyperparameters
  piw_structure = c("genealogy", "institution"),
  piw_skill_weight = 0.5,
  piw_transform = c("gaussian", "softmax", "rank", "linear"),
  piw_metric = c("rmse", "correlation", "bias", "nse"),
  piw_temperature = NULL,
  piw_robust_scale = c("sd", "mad"),
  min_weight = 1e-3,
  verbose = FALSE
) {
  piw_structure <- match.arg(piw_structure)
  piw_transform <- match.arg(piw_transform)
  piw_metric <- match.arg(piw_metric)
  piw_robust_scale <- match.arg(piw_robust_scale)

  # --- 1) Equal weights (uniform structural, no skill) ---
  w_equal <- compute_gcm_weights_bma(
    gcm_data = gcm_data,
    skill_data = skill_data,
    kuma_table = kuma_table,
    model_col = model_col,
    scenario_col = scenario_col,
    skill_col = skill_col,
    prior_method = "uniform",
    skill_metric = piw_metric,
    likelihood_transform = piw_transform,
    skill_weight = 0,          # structure-only => uniform
    min_weight = min_weight,
    temperature = piw_temperature,
    robust_scale = piw_robust_scale,
    include_missing_skill = FALSE,
    verbose = verbose,
    skill_weight_grid = NULL
  )

  # --- 2) Institution weights (structure-only) ---
  w_inst <- compute_gcm_weights_bma(
    gcm_data = gcm_data,
    skill_data = skill_data,
    kuma_table = kuma_table,
    model_col = model_col,
    scenario_col = scenario_col,
    skill_col = skill_col,
    prior_method = "institution",
    institution_col = institution_col,
    skill_metric = piw_metric,
    likelihood_transform = piw_transform,
    skill_weight = 0,          # structure-only
    min_weight = min_weight,
    temperature = piw_temperature,
    robust_scale = piw_robust_scale,
    include_missing_skill = FALSE,
    verbose = verbose,
    skill_weight_grid = NULL
  )

  # --- 3) Genealogy weights (structure-only) ---
  w_gen <- compute_gcm_weights_bma(
    gcm_data = gcm_data,
    skill_data = skill_data,
    kuma_table = kuma_table,
    model_col = model_col,
    scenario_col = scenario_col,
    skill_col = skill_col,
    prior_method = "genealogy",
    skill_metric = piw_metric,
    likelihood_transform = piw_transform,
    skill_weight = 0,          # structure-only
    min_weight = min_weight,
    temperature = piw_temperature,
    robust_scale = piw_robust_scale,
    include_missing_skill = FALSE,
    verbose = verbose,
    skill_weight_grid = NULL
  )

  # --- 4) PIW (structure + performance) ---
  w_piw <- compute_gcm_weights_bma(
    gcm_data = gcm_data,
    skill_data = skill_data,
    kuma_table = kuma_table,
    model_col = model_col,
    scenario_col = scenario_col,
    skill_col = skill_col,
    prior_method = piw_structure,
    institution_col = institution_col,
    skill_metric = piw_metric,
    likelihood_transform = piw_transform,
    skill_weight = piw_skill_weight,
    min_weight = min_weight,
    temperature = piw_temperature,
    robust_scale = piw_robust_scale,
    include_missing_skill = FALSE,
    verbose = verbose,
    skill_weight_grid = seq(0, 1, by = 0.25)
  )

  # Stats per scenario (weighted distribution over rows in gcm_data)
  s_equal <- compute_weighted_ensemble_stats(gcm_data, w_equal, model_col, scenario_col, target_col, weight_col = "w_bma")
  s_inst  <- compute_weighted_ensemble_stats(gcm_data, w_inst,  model_col, scenario_col, target_col, weight_col = "w_bma")
  s_gen   <- compute_weighted_ensemble_stats(gcm_data, w_gen,   model_col, scenario_col, target_col, weight_col = "w_bma")
  s_piw   <- compute_weighted_ensemble_stats(gcm_data, w_piw,   model_col, scenario_col, target_col, weight_col = "w_bma")

  # Attach method label
  s_equal$method <- "equal"
  s_inst$method  <- "institution"
  s_gen$method   <- "genealogy"
  s_piw$method   <- paste0("piw_", piw_structure)

  stats_all <- rbind(s_equal, s_inst, s_gen, s_piw)
  stats_all <- stats_all[order(stats_all$scenario, stats_all$method), , drop = FALSE]

  list(
    weights = list(equal = w_equal, institution = w_inst, genealogy = w_gen, piw = w_piw),
    stats = stats_all
  )
}

# =============================================================================
# 3) Run synthetic example + compare findings
# =============================================================================


ex <- make_synthetic_gcm_example(seed = 7)

out <- compute_all_weighting_options(
  gcm_data = ex$gcm_future,
  skill_data = ex$skill_df,
  kuma_table = ex$kuma_table,
  piw_structure = "genealogy",
  piw_skill_weight = 0.6,      # emphasize performance moderately
  piw_transform = "gaussian",  # monotone to best-score
  piw_metric = "rmse",
  piw_temperature = NULL,
  piw_robust_scale = "mad",
  verbose = TRUE
)

# --- Summary table: weighted mean/SD/quantiles by scenario and method
print(out$stats)

# --- Compare weight concentration for each method (per scenario)
summarize_weight_concentration <- function(weights_df_list) {
  methods <- names(weights_df_list)
  res <- list()
  k <- 1
  for (m in methods) {
    wdf <- weights_df_list[[m]]
    for (sc in sort(unique(wdf$scenario))) {
      w <- wdf$w_bma[wdf$scenario == sc]
      w <- w[is.finite(w)]
      if (length(w) == 0) next
      n_eff <- 1 / sum(w^2)
      top <- max(w)
      res[[k]] <- data.frame(
        method = m,
        scenario = sc,
        n_models = length(w),
        n_eff = n_eff,
        top_weight = top,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  do.call(rbind, res)
}

conc <- summarize_weight_concentration(out$weights)
print(conc)

# --- Optional: validate sums (from utils)
validate_weight_sums(out$weights$piw, scenario_col = "scenario", weight_col = "w_bma")
