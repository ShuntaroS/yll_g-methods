# Core estimation engine.
#
# Performs one YLL estimation under arbitrary `intervention_reference` /
# `intervention_exposed` specifications and an optional `target_population`.
# This is the function that all user-facing wrappers (`estimate_yll_gformula`,
# `*_intervention`, `*_binary_stochastic`, the ATE/ATT/ATC wrappers, ...)
# eventually call. Doing one bootstrap iteration also calls this function.
#
# The pipeline is:
#
#   1. Validate inputs and discretise follow-up time.
#   2. Fit the discrete-time pooled logistic hazard model on the observed
#      data (full sample — see the note in `intervention.R` about why we do
#      not subset by estimand here).
#   3. Build a counterfactual prediction grid (one row per subject per
#      integer age in [age_start, age_end]).
#   4. Resolve the two intervention specifications into per-row P(A* = exp).
#   5. Predict per-row hazards under both arms, turn them into per-subject
#      survival curves, then average across subjects within the chosen
#      target population to get the marginal counterfactual curves.
#   6. For every requested starting age a_start, condition on survival to
#      a_start and integrate the residual life expectancies under each arm;
#      YLL is the difference (LE_ref - LE_exp).
#
# Returns a tibble with one row per starting age, plus the per-arm marginal
# survival curve as an attribute (used by the bootstrap collector and the
# plotting helpers).
#' @noRd
estimate_yll_gformula_engine_single <- function(
    data,
    id_var,
    time_var,
    event_var,
    exposure_var,
    reference_level,
    exposed_level,
    age_at_entry_var,
    age_start,
    age_end,
    age_interval,
    confounders_baseline = NULL,
    intervention_reference,
    intervention_exposed,
    target_population = NULL,
    cut_by = 1L,
    integration = c("left_rectangle", "trapezoidal")
) {
  integration <- match.arg(integration)

  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )

  cut_points <- yll_make_cut_points(data, time_var = time_var, by = cut_by)

  # Person-period long form, then derive `age_temp` (= attained age within
  # the interval) and `expo` (a stable name for the exposure) so the rest of
  # the pipeline doesn't have to know the user's column names.
  df_long <- yll_survsplit(
    data       = data,
    time_var   = time_var,
    event_var  = event_var,
    cut_points = cut_points
  ) |>
    mutate(
      age_temp = .data[[age_at_entry_var]] + tstart,
      expo     = .data[[exposure_var]]
    )

  knot_p <- yll_knot_quantiles(df_long$age_temp)

  fit <- yll_fit_discrete_hazard_glm(
    df_long              = df_long,
    event_var            = event_var,
    confounders_baseline = confounders_baseline,
    knot_p               = knot_p
  )

  # Counterfactual data: only build rows for ages in [age_start, age_end].
  #
  # Mathematically we could expand to the full lifetime, but the conditional
  # survival S(t)/S(a_start) cancels every hazard contribution at ages
  # smaller than a_start, so restricting to this window gives identical
  # results while (i) reducing compute roughly proportionally and (ii)
  # avoiding extrapolation of the hazard model into ages outside the
  # observed support.
  cf_base <- yll_expand_counterfactual_data(
    data           = data,
    id_var         = id_var,
    exposure_var   = exposure_var,
    age_temp_start = as.integer(age_start),
    age_temp_end   = as.integer(age_end)
  )

  # Per-subject P(A* = exposed) under each intervention arm. These are vectors
  # of length nrow(cf_base) and will be used to mix the two per-arm survival
  # curves into the marginal counterfactual curves.
  p_exposed_reference <- yll_resolve_exposed_probability(
    cf_base,
    intervention_reference,
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  p_exposed_exposed <- yll_resolve_exposed_probability(
    cf_base,
    intervention_exposed,
    reference_level = reference_level,
    exposed_level = exposed_level
  )

  # Predict hazards under "everyone reference" and "everyone exposed" copies
  # of the counterfactual data, then turn them into per-subject survival
  # curves via cumulative product.
  cf0 <- cf_base
  cf1 <- cf_base
  cf0$expo <- reference_level
  cf1$expo <- exposed_level
  cf_base$hazard0 <- yll_predict_hazard(fit, cf0)
  cf_base$hazard1 <- yll_predict_hazard(fit, cf1)
  cf_base <- yll_compute_individual_survival_curves(cf_base, id_var = id_var)

  # Apply the target-population mask *after* counterfactual curves are built.
  # This way the hazard model is always fit on the full sample, and only the
  # averaging step changes when the user selects a different estimand.
  pop_idx <- yll_resolve_population_index(cf_base, target_population)
  cf_pop <- cf_base[pop_idx, , drop = FALSE]
  p_exposed_reference <- p_exposed_reference[pop_idx]
  p_exposed_exposed <- p_exposed_exposed[pop_idx]

  s0 <- yll_mean_survival_under_intervention(cf_pop, p_exposed_reference, surv_name = "surv0")
  s1 <- yll_mean_survival_under_intervention(cf_pop, p_exposed_exposed, surv_name = "surv1")

  # `g_surv` holds the two marginal counterfactual survival curves on the
  # same age grid; we keep it as one tibble for easy downstream use.
  g_surv <- inner_join(s0, s1, by = "age_temp")

  age_list <- seq(from = age_start, to = age_end, by = age_interval)

  # For each requested starting age: condition the marginal curves on
  # survival to that age, integrate to get residual LE under each arm, and
  # report the difference as YLL.
  yll_table <- map_dfr(age_list, function(a_start) {
    df_cond <- yll_conditional_survival_from_age(g_surv, a_start)
    if (is.null(df_cond)) {
      return(tibble(age_start = a_start, yll = NA_real_, le_m0 = NA_real_, le_m1 = NA_real_))
    }

    le_m0 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g0, integration = integration)
    le_m1 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g1, integration = integration)

    tibble(
      age_start = a_start,
      yll  = le_m0 - le_m1,
      le_m0 = le_m0,
      le_m1 = le_m1
    )
  })

  # Stash the marginal curves on the returned tibble so the bootstrap
  # collector and plotting helpers can pick them up without recomputing.
  attr(yll_table, "marginal_curves") <- as_tibble(g_surv)
  yll_table
}

# Backward-compatible wrapper around the engine for the named-estimand
# (ATE / ATT / ATC) interface.
#
# It resolves binary exposure levels, then turns the estimand into a
# `target_population` rule and forwards everything to the engine. The
# `intervention_reference` / `intervention_exposed` arguments are filled in
# with the deterministic "everyone reference" vs "everyone exposed" pair —
# which is exactly the contrast the named estimands require.
#' @noRd
estimate_yll_gformula_single <- function(
    data,
    id_var,
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start,
    age_end,
    age_interval,
    confounders_baseline = NULL,
    estimand = c("ATT", "ATC", "ATE"),
    cut_by = 1L,
    integration = c("left_rectangle", "trapezoidal")
) {
  estimand <- match.arg(estimand)
  integration <- match.arg(integration)

  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level

  estimate_yll_gformula_engine_single(
    data = data,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    exposure_var = exposure_var,
    reference_level = reference_level,
    exposed_level = exposed_level,
    age_at_entry_var = age_at_entry_var,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    confounders_baseline = confounders_baseline,
    intervention_reference = reference_level,
    intervention_exposed = exposed_level,
    target_population = yll_make_population_rule_from_estimand(
      estimand = estimand,
      reference_level = reference_level,
      exposed_level = exposed_level
    ),
    cut_by = cut_by,
    integration = integration
  )
}
