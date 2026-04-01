load_legacy_functions <- function() {
  legacy_code <- system2("git", c("show", "HEAD:yll_g_formula.R"), stdout = TRUE)
  legacy_env <- new.env(parent = globalenv())
  eval(parse(text = legacy_code), envir = legacy_env)
  legacy_env
}

simulate_fixed_entry_data <- function(seed = 20260401, n = 3000) {
  set.seed(seed)
  
  max_follow <- 8
  age0 <- rep(50, n)
  expo <- factor(rbinom(n, 1, 0.45), levels = c(0, 1), labels = c("Never", "Current/Ever"))
  c_age <- rnorm(n)
  
  alpha <- -4.5
  beta_a <- 0.5
  beta_c <- 0.3
  beta_age <- 0.18
  
  period <- rep(max_follow, n)
  dflag <- integer(n)
  
  for (i in seq_len(n)) {
    for (j in 0:(max_follow - 1)) {
      hazard_ij <- plogis(alpha + beta_a * (expo[i] == "Current/Ever") + beta_c * c_age[i] + beta_age * j)
      if (runif(1) < hazard_ij) {
        period[i] <- j + 1
        dflag[i] <- 1L
        break
      }
    }
  }
  
  data.frame(
    id = seq_len(n),
    period = period,
    dflag = dflag,
    smoke_binary = expo,
    age = age0,
    z = c_age
  )
}

as_legacy_point <- function(res) {
  res$detailed_results |>
    rename(
      age_start = starting_age,
      le_m0 = le_reference,
      le_m1 = le_exposed
    ) |>
    select(age_start, yll, le_m0, le_m1)
}

as_legacy_summary <- function(res) {
  out <- res$detailed_results
  
  missing_cols <- setdiff(
    c(
      "yll_se_normal",
      "le_reference_se_normal",
      "le_exposed_se_normal",
      "ci_low_percentile",
      "ci_high_percentile",
      "le_reference_ci_low_percentile",
      "le_reference_ci_high_percentile",
      "le_exposed_ci_low_percentile",
      "le_exposed_ci_high_percentile",
      "ci_low_normal",
      "ci_high_normal",
      "le_reference_ci_low_normal",
      "le_reference_ci_high_normal",
      "le_exposed_ci_low_normal",
      "le_exposed_ci_high_normal"
    ),
    names(out)
  )
  
  for (col in missing_cols) {
    out[[col]] <- NA_real_
  }
  
  out |>
    rename(
      age_start = starting_age,
      le_m0 = le_reference,
      le_m1 = le_exposed,
      se_yll_norm = yll_se_normal,
      se_le_m0_norm = le_reference_se_normal,
      se_le_m1_norm = le_exposed_se_normal,
      yll_lwr_perc = ci_low_percentile,
      yll_upr_perc = ci_high_percentile,
      le_m0_lwr_perc = le_reference_ci_low_percentile,
      le_m0_upr_perc = le_reference_ci_high_percentile,
      le_m1_lwr_perc = le_exposed_ci_low_percentile,
      le_m1_upr_perc = le_exposed_ci_high_percentile,
      yll_lwr_norm = ci_low_normal,
      yll_upr_norm = ci_high_normal,
      le_m0_lwr_norm = le_reference_ci_low_normal,
      le_m0_upr_norm = le_reference_ci_high_normal,
      le_m1_lwr_norm = le_exposed_ci_low_normal,
      le_m1_upr_norm = le_exposed_ci_high_normal
    )
}

assert_same_result <- function(x, y, label) {
  y_summary_legacy <- as_legacy_summary(y)
  y_summary_legacy <- y_summary_legacy[, names(x$summary), drop = FALSE]
  
  stopifnot(isTRUE(all.equal(x$point, as_legacy_point(y), tolerance = 0, check.attributes = FALSE)))
  stopifnot(isTRUE(all.equal(x$summary, y_summary_legacy, tolerance = 0, check.attributes = FALSE)))
  cat(label, "matched exactly.\n")
}

assert_same_current_result <- function(x, y, label) {
  stopifnot(isTRUE(all.equal(x$detailed_results, y$detailed_results, tolerance = 0, check.attributes = FALSE)))
  stopifnot(isTRUE(all.equal(x$summary, y$summary, tolerance = 0, check.attributes = FALSE)))
  cat(label, "matched exactly.\n")
}

legacy <- load_legacy_functions()
source("yll_g_formula.R")

sim_data <- simulate_fixed_entry_data()

common_args <- list(
  data = sim_data,
  B = 5,
  seed = 11,
  conf_level = 0.95,
  method = "normal",
  show_progress = FALSE,
  id_var = "id",
  time_var = "period",
  event_var = "dflag",
  exposure_var = "smoke_binary",
  reference_level = "Never",
  exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 58,
  age_interval = 4,
  confounders_baseline = "z",
  use_future = FALSE
)

intervention_args <- common_args
intervention_args$reference_level <- NULL
intervention_args$exposed_level <- NULL

legacy_ate <- do.call(legacy$estimate_yll_gformula, c(common_args, list(estimand = "ATE")))
current_ate <- do.call(estimate_yll_gformula, c(common_args, list(estimand = "ATE")))
assert_same_result(legacy_ate, current_ate, "Refactored estimate_yll_gformula(ATE)")

current_att <- do.call(estimate_yll_gformula, c(common_args, list(estimand = "ATT")))
wrapper_att <- do.call(estimate_yll_gformula_att, common_args)
assert_same_current_result(current_att, wrapper_att, "estimate_yll_gformula_att()")

current_atc <- do.call(estimate_yll_gformula, c(common_args, list(estimand = "ATC")))
wrapper_atc <- do.call(estimate_yll_gformula_atc, common_args)
assert_same_current_result(current_atc, wrapper_atc, "estimate_yll_gformula_atc()")

wrapper_ate <- do.call(estimate_yll_gformula_ate, common_args)
assert_same_current_result(current_ate, wrapper_ate, "estimate_yll_gformula_ate()")

generic_ate <- do.call(
  estimate_yll_gformula_intervention,
  c(
    intervention_args,
    list(
      intervention_reference = "Never",
      intervention_exposed = "Current/Ever",
      target_population = NULL
    )
  )
)
assert_same_current_result(current_ate, generic_ate, "estimate_yll_gformula_intervention() ATE")

generic_att <- do.call(
  estimate_yll_gformula_intervention,
  c(
    intervention_args,
    list(
      intervention_reference = "Never",
      intervention_exposed = "Current/Ever",
      target_population = function(data) data$expo_original == "Current/Ever"
    )
  )
)
assert_same_current_result(current_att, generic_att, "estimate_yll_gformula_intervention() ATT-style population")

stochastic_shift <- do.call(
  estimate_yll_gformula_intervention,
  c(
    modifyList(intervention_args, list(B = 0, method = "normal")),
    list(
      intervention_reference = function(data) data$expo_original,
      intervention_exposed = function(data) ifelse(data$expo_original == "Current/Ever", 0.7, 0),
      target_population = NULL
    )
  )
)

stopifnot(nrow(stochastic_shift$detailed_results) > 0)
stopifnot(all(is.finite(stochastic_shift$detailed_results$yll)))
cat("estimate_yll_gformula_intervention() stochastic shift ran successfully.\n")

stochastic_increase_among_unexposed <- do.call(
  estimate_yll_gformula_intervention,
  c(
    modifyList(intervention_args, list(B = 0, method = "normal")),
    list(
      intervention_reference = function(data) data$expo_original,
      intervention_exposed = yll_make_binary_stochastic_intervention(
        prob_exposed_if_unexposed = 0.2,
        prob_exposed_if_exposed = 1,
        reference_level = "Never",
        exposed_level = "Current/Ever"
      ),
      target_population = NULL
    )
  )
)

stopifnot(nrow(stochastic_increase_among_unexposed$detailed_results) > 0)
stopifnot(all(is.finite(stochastic_increase_among_unexposed$detailed_results$yll)))
cat("estimate_yll_gformula_intervention() exposure increase among unexposed ran successfully.\n")

wrapper_stochastic_shift <- do.call(
  estimate_yll_gformula_binary_stochastic_vs_natural,
  c(
    modifyList(common_args, list(B = 0, method = "normal")),
    list(
      prob_exposed_if_unexposed = 0,
      prob_exposed_if_exposed = 0.7
    )
  )
)
assert_same_current_result(
  stochastic_shift,
  wrapper_stochastic_shift,
  "estimate_yll_gformula_binary_stochastic_vs_natural() exposure decrease among exposed"
)

wrapper_stochastic_increase_among_unexposed <- do.call(
  estimate_yll_gformula_binary_stochastic_vs_natural,
  c(
    modifyList(common_args, list(B = 0, method = "normal")),
    list(
      prob_exposed_if_unexposed = 0.2,
      prob_exposed_if_exposed = 1
    )
  )
)
assert_same_current_result(
  stochastic_increase_among_unexposed,
  wrapper_stochastic_increase_among_unexposed,
  "estimate_yll_gformula_binary_stochastic_vs_natural() exposure increase among unexposed"
)

cat("\nAll equivalence checks passed.\n")
