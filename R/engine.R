# 共通エンジン本体: 任意の介入と対象集団の下で 1 回分の YLL 推定を行う。
# Core engine: performs one YLL estimation under arbitrary interventions and target population.
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

  # 反実仮想データ: 年齢 [age_start, age_end] のみを作る。
  # 条件付き生存 S(t)/S(a_start) は a_start より小さい年齢の hazard をキャンセルするため、
  # この範囲だけ展開すれば結果は同じで、計算量が減りモデル外挿のリスクも下がる。
  cf_base <- yll_expand_counterfactual_data(
    data           = data,
    id_var         = id_var,
    exposure_var   = exposure_var,
    age_temp_start = as.integer(age_start),
    age_temp_end   = as.integer(age_end)
  )

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

  cf0 <- cf_base
  cf1 <- cf_base
  cf0$expo <- reference_level
  cf1$expo <- exposed_level
  cf_base$hazard0 <- yll_predict_hazard(fit, cf0)
  cf_base$hazard1 <- yll_predict_hazard(fit, cf1)
  cf_base <- yll_compute_individual_survival_curves(cf_base, id_var = id_var)

  pop_idx <- yll_resolve_population_index(cf_base, target_population)
  cf_pop <- cf_base[pop_idx, , drop = FALSE]
  p_exposed_reference <- p_exposed_reference[pop_idx]
  p_exposed_exposed <- p_exposed_exposed[pop_idx]

  s0 <- yll_mean_survival_under_intervention(cf_pop, p_exposed_reference, surv_name = "surv0")
  s1 <- yll_mean_survival_under_intervention(cf_pop, p_exposed_exposed, surv_name = "surv1")

  g_surv <- inner_join(s0, s1, by = "age_temp")

  age_list <- seq(from = age_start, to = age_end, by = age_interval)

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

  attr(yll_table, "marginal_curves") <- as_tibble(g_surv)
  yll_table
}

# 既存仕様ラッパー: ATT/ATC/ATE 指定を共通エンジンへ橋渡しする。
# Backward-compatible wrapper: maps ATT/ATC/ATE inputs to the shared engine.
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
