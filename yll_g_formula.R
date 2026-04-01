library(tidyverse)
library(Epi)      # Ns()
library(survival) # Surv(), survSplit()
library(future)
library(future.apply)
library(progressr)


# カット点を作る: 離散時間化のための分割点を追跡時間から作成する。
# Build cut points: creates split points for discretizing follow-up time.
yll_make_cut_points <- function(data, time_var, by = 1L) {
  # survSplit の cut は 0 を含めないベクトルを渡す想定なので，ここでは全体を返す
  max_t <- ceiling(max(data[[time_var]], na.rm = TRUE))
  seq.int(from = 0L, to = max_t, by = by)
}

# 生存データをロング形式へ変換する: `survSplit()` で時間区間ごとのデータに展開する。
# Split survival data: expands person-level data into interval-level rows with `survSplit()`.
yll_survsplit <- function(data, time_var, event_var, cut_points) {
  # Surv(time, event) ~ . を survSplit に渡す
  surv_formula <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ."))
  survSplit(
    formula  = surv_formula,
    data     = data,
    cut      = cut_points[-1],
    episode  = "tgroup"
  )
}

# ノット候補を作る: 年齢変数の分位点から spline のノットを作成する。
# Build spline knots: uses quantiles of age to define spline knots.
yll_knot_quantiles <- function(x, probs = c(0, 0.25, 0.50, 0.75, 1)) {
  as.numeric(quantile(x, probs = probs, na.rm = TRUE))
}

# 交絡因子の右辺項を作る: モデル式に入れる共変量部分を文字列で返す。
# Build confounder RHS: returns the covariate part of the model formula.
yll_rhs_confounders <- function(confounders_baseline) {
  if (is.null(confounders_baseline) || length(confounders_baseline) == 0) {
    "1"
  } else {
    paste(confounders_baseline, collapse = " + ")
  }
}

# 必須列を確認する: 指定された変数名が `data` に存在するか検査する。
# Check required columns: verifies that requested variables exist in `data`.
yll_check_required_columns <- function(data, vars) {
  missing_vars <- setdiff(unique(vars), names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Missing columns in `data`: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
}

# 二値曝露のレベルを決める: 参照群と曝露群の値を確定する。
# Resolve binary exposure levels: determines reference and exposed values.
yll_resolve_binary_levels <- function(x, reference_level = NULL, exposed_level = NULL) {
  observed <- unique(x[!is.na(x)])
  
  if (length(observed) != 2) {
    stop("Currently only binary exposures are supported.", call. = FALSE)
  }
  
  inferred_levels <- if (is.factor(x)) levels(x) else observed
  inferred_levels <- inferred_levels[inferred_levels %in% observed]
  
  if (is.null(reference_level)) {
    reference_level <- inferred_levels[[1]]
  }
  if (is.null(exposed_level)) {
    exposed_level <- inferred_levels[[2]]
  }
  
  if (!(reference_level %in% observed)) {
    stop("`reference_level` was not found in the observed exposure values.", call. = FALSE)
  }
  if (!(exposed_level %in% observed)) {
    stop("`exposed_level` was not found in the observed exposure values.", call. = FALSE)
  }
  if (identical(reference_level, exposed_level)) {
    stop("`reference_level` and `exposed_level` must be different.", call. = FALSE)
  }
  
  list(reference_level = reference_level, exposed_level = exposed_level)
}

# 離散時間ハザードモデルを当てる: 年齢 spline と曝露の交互作用を含む GLM を推定する。
# Fit discrete-time hazard model: estimates a GLM with age splines and exposure interaction.
yll_fit_discrete_hazard_glm <- function(
    df_long,
    event_var,
    confounders_baseline = NULL,
    knot_p
) {
  rhs <- yll_rhs_confounders(confounders_baseline)
  
  # 交互作用込み
  form <- as.formula(
    paste0(
      event_var,
      " ~ ",
      "expo",
      " + Epi::Ns(age_temp, knots = knot_p) + ",
      "expo * Epi::Ns(age_temp, knots = knot_p) + ",
      rhs
    )
  )
  
  glm(formula = form, family = binomial(), data = df_long)
}

# 反実仮想データを展開する: 各個体について年齢ごとの予測用データを作る。
# Expand counterfactual data: creates age-specific prediction rows for each individual.
yll_expand_counterfactual_data <- function(
    data,
    id_var,
    exposure_var,
    age_end_inclusive
) {
  # 各個体について age_temp = 0..age_end_inclusive を作成
  # uncount を使うため tidyr 前提
  out <- data |>
    uncount(weights = age_end_inclusive + 1L, .remove = FALSE) |>
    group_by(.data[[id_var]]) |>
    mutate(age_temp = row_number() - 1L) |>
    ungroup() |>
    rename(expo_original = all_of(exposure_var))
  
  # モデル用 exposure 列は別名で持つ
  out[["expo"]] <- out$expo_original
  out
}

# ハザードを予測する: 学習済みモデルから各行の死亡確率を返す。
# Predict hazards: returns row-wise event probabilities from the fitted model.
yll_predict_hazard <- function(fit, newdata) {
  as.numeric(predict(fit, newdata = newdata, type = "response"))
}

# 対象集団ルールを作る: ATT/ATC/ATE を logical な選択規則へ変換する。
# Build population rule: converts ATT/ATC/ATE into a logical selection rule.
yll_make_population_rule_from_estimand <- function(estimand, reference_level, exposed_level) {
  if (identical(estimand, "ATE")) {
    return(NULL)
  }
  
  if (identical(estimand, "ATT")) {
    return(function(data) data$expo_original == exposed_level)
  }
  
  if (identical(estimand, "ATC")) {
    return(function(data) data$expo_original == reference_level)
  }
  
  stop("estimand must be one of 'ATT', 'ATC', or 'ATE'.", call. = FALSE)
}

# 対象集団を解決する: 平均化の対象となる個体を logical ベクトルで返す。
# Resolve target population: returns a logical index for the population to average over.
yll_resolve_population_index <- function(data, target_population = NULL) {
  if (is.null(target_population)) {
    return(rep(TRUE, nrow(data)))
  }
  
  idx <- if (is.function(target_population)) target_population(data) else target_population
  
  if (!is.logical(idx)) {
    stop("`target_population` must evaluate to a logical vector.", call. = FALSE)
  }
  if (length(idx) != nrow(data)) {
    stop("`target_population` must have length equal to `nrow(data)`.", call. = FALSE)
  }
  if (anyNA(idx)) {
    stop("`target_population` must not contain `NA` values.", call. = FALSE)
  }
  
  idx
}

# 介入確率を解決する: 各個体の `A*=1` の確率へ介入指定を変換する。
# Resolve intervention probabilities: converts an intervention specification into subject-specific probabilities of exposure.
yll_resolve_exposed_probability <- function(data, intervention, reference_level, exposed_level) {
  assigned <- if (is.function(intervention)) intervention(data) else intervention
  
  if (length(assigned) == 1L) {
    assigned <- rep(assigned, nrow(data))
  }
  if (length(assigned) != nrow(data)) {
    stop("Intervention must return a scalar or a vector of length `nrow(data)`.", call. = FALSE)
  }
  
  if (is.numeric(assigned)) {
    if (anyNA(assigned) || any(assigned < 0 | assigned > 1)) {
      stop("Numeric interventions must be probabilities between 0 and 1.", call. = FALSE)
    }
    return(as.numeric(assigned))
  }
  
  if (is.logical(assigned)) {
    if (anyNA(assigned)) {
      stop("Logical interventions must not contain `NA` values.", call. = FALSE)
    }
    return(as.numeric(assigned))
  }
  
  assigned_chr <- as.character(assigned)
  reference_chr <- as.character(reference_level)
  exposed_chr <- as.character(exposed_level)
  
  if (!all(assigned_chr %in% c(reference_chr, exposed_chr))) {
    stop(
      "Intervention values must be either probabilities in [0, 1] or exposure values matching the observed binary levels.",
      call. = FALSE
    )
  }
  
  as.numeric(assigned_chr == exposed_chr)
}

# 二値曝露の stochastic 介入を作る: 現在の曝露状態ごとに `A*=1` の確率を指定する。
# Build binary stochastic intervention: specifies `P(A*=1)` separately for currently unexposed and exposed individuals.
yll_make_binary_stochastic_intervention <- function(
    prob_exposed_if_unexposed,
    prob_exposed_if_exposed,
    reference_level,
    exposed_level
) {
  probs <- c(prob_exposed_if_unexposed, prob_exposed_if_exposed)
  if (anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop("Intervention probabilities must be between 0 and 1.", call. = FALSE)
  }
  
  reference_chr <- as.character(reference_level)
  exposed_chr <- as.character(exposed_level)
  
  function(data) {
    expo_chr <- as.character(data$expo_original)
    
    if (!all(expo_chr %in% c(reference_chr, exposed_chr))) {
      stop("`expo_original` contains values outside the specified binary exposure levels.", call. = FALSE)
    }
    
    ifelse(
      expo_chr == reference_chr,
      prob_exposed_if_unexposed,
      prob_exposed_if_exposed
    )
  }
}

# 個体別生存を作る: 参照曝露と曝露群の下での生存曲線を個体ごとに計算する。
# Build individual survival curves: computes per-subject survival under the reference and exposed levels.
yll_compute_individual_survival_curves <- function(df, id_var, hazard0_var = "hazard0", hazard1_var = "hazard1") {
  df |>
    arrange(.data[[id_var]], age_temp) |>
    group_by(.data[[id_var]]) |>
    mutate(
      .surv0_after_interval = cumprod(1 - .data[[hazard0_var]]),
      .surv1_after_interval = cumprod(1 - .data[[hazard1_var]]),
      surv0 = dplyr::lag(.surv0_after_interval, default = 1),
      surv1 = dplyr::lag(.surv1_after_interval, default = 1)
    ) |>
    ungroup()
}

# 介入下の平均生存を作る: 個体別の曝露確率で2本の生存曲線を混合して平均する。
# Build intervention survival: mixes the two survival curves using subject-specific exposure probabilities.
yll_mean_survival_under_intervention <- function(df, p_exposed, surv_name = "surv") {
  if (length(p_exposed) != nrow(df)) {
    stop("`p_exposed` must have length equal to `nrow(df)`.", call. = FALSE)
  }
  
  df |>
    mutate(
      .p_exposed = p_exposed,
      .surv_mix = (1 - .p_exposed) * surv0 + .p_exposed * surv1
    ) |>
    group_by(age_temp) |>
    summarise(!!surv_name := mean(.surv_mix), .groups = "drop")
}

# 年齢別平均生存を作る: 個体ごとの生存を計算し、各年齢で平均する。
# Average survival by age: computes individual survival and averages it at each age.
yll_mean_survival_by_age <- function(df, id_var, hazard_var, surv_name = "surv") {
  # Survival is indexed at the start of each age interval.
  df |>
    arrange(.data[[id_var]], age_temp) |>
    group_by(.data[[id_var]]) |>
    mutate(
      .surv_after_interval = cumprod(1 - .data[[hazard_var]]),
      .surv_at_age = dplyr::lag(.surv_after_interval, default = 1)
    ) |>
    ungroup() |>
    group_by(age_temp) |>
    summarise(!!surv_name := mean(.surv_at_age), .groups = "drop")
}

# 条件付き生存を作る: 開始年齢で標準化した生存曲線に変換する。
# Build conditional survival: rescales survival curves to start at a chosen age.
yll_conditional_survival_from_age <- function(g_surv_data, a_start) {
  s0 <- g_surv_data[["surv0"]][g_surv_data$age_temp == a_start]
  s1 <- g_surv_data[["surv1"]][g_surv_data$age_temp == a_start]
  
  if (length(s0) == 0 || length(s1) == 0 || is.na(s0) || is.na(s1) || s0 <= 0 || s1 <= 0) {
    return(NULL)
  }
  
  g_surv_data |>
    filter(age_temp >= a_start) |>
    mutate(
      surv_cond_g0 = .data[["surv0"]] / s0,
      surv_cond_g1 = .data[["surv1"]] / s1
    )
}

# 余命を積分する: 条件付き生存曲線の面積から LE を計算する。
# Integrate life expectancy: computes LE as the area under the conditional survival curve.
yll_integrate_le <- function(age, surv) {
  if (length(age) <= 1) {
    return(0)
  }
  
  ord <- order(age)
  age <- age[ord]
  surv <- surv[ord]
  
  sum(diff(age) * head(surv, -1))
}

# 共通エンジン本体: 任意の介入と対象集団の下で 1 回分の YLL 推定を行う。
# Core engine: performs one YLL estimation under arbitrary interventions and target population.
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
    cut_by = 1L
) {
  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  
  # survSplit 用のカット
  cut_points <- yll_make_cut_points(data, time_var = time_var, by = cut_by)
  
  # 離散時間ロング化
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
  
  # 反実仮想データ（年齢 0..age_end を作る）
  age_end_inclusive <- as.integer(age_end)
  cf_base <- yll_expand_counterfactual_data(
    data              = data,
    id_var            = id_var,
    exposure_var      = exposure_var,
    age_end_inclusive = age_end_inclusive
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
  
  # 生存曲線（年齢ごと平均）
  s0 <- yll_mean_survival_under_intervention(cf_pop, p_exposed_reference, surv_name = "surv0")
  s1 <- yll_mean_survival_under_intervention(cf_pop, p_exposed_exposed, surv_name = "surv1")
  
  g_surv <- inner_join(s0, s1, by = "age_temp")
  
  age_list <- seq(from = age_start, to = age_end, by = age_interval)
  
  map_dfr(age_list, function(a_start) {
    df_cond <- yll_conditional_survival_from_age(g_surv, a_start)
    if (is.null(df_cond)) {
      return(tibble(age_start = a_start, yll = NA_real_, le_m0 = NA_real_, le_m1 = NA_real_))
    }
    
    le_m0 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g0)
    le_m1 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g1)
    
    tibble(
      age_start = a_start,
      yll  = le_m0 - le_m1,
      le_m0 = le_m0,
      le_m1 = le_m1
    )
  })
}

# 既存仕様ラッパー: ATT/ATC/ATE 指定を共通エンジンへ橋渡しする。
# Backward-compatible wrapper: maps ATT/ATC/ATE inputs to the shared engine.
etimate_yll_gformula_single <- function(
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
    cut_by = 1L
) {
  estimand <- match.arg(estimand)
  
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
    cut_by = cut_by
  )
}

# ------------------------------------------------------------------------------
# Bootstrap
# ------------------------------------------------------------------------------

# bootstrap 用 ID を標本抽出する: 個体 ID を復元抽出する。
# Sample bootstrap IDs: resamples individual IDs with replacement.
yll_bootstrap_ids <- function(ids) {
  n <- length(ids)
  sample(ids, size = n, replace = TRUE)
}

# bootstrap データを作る: 重複抽出された個体に新しい ID を振る。
# Build bootstrap data: assigns fresh IDs to duplicated sampled individuals.
yll_make_boot_data <- function(data, id_var, sampled_ids) {
  # 同一 id が複数回サンプルされた場合に備えて一意な boot_id を振る
  id_map <- tibble(
    !!id_var := sampled_ids,
    boot_id  = paste0(sampled_ids, "_rep", ave(sampled_ids, sampled_ids, FUN = seq_along))
  )
  
  left_join(id_map, data, by = id_var) |>
    mutate(!!id_var := boot_id) |>
    select(-boot_id)
}

# パーセンタイル CI を作る: bootstrap 分布から信頼区間を計算する。
# Compute percentile CIs: derives confidence intervals from bootstrap quantiles.
yll_ci_percentile <- function(boot_df, conf_level) {
  alpha <- (1 - conf_level) / 2
  boot_df |>
    group_by(age_start) |>
    summarise(
      yll_lwr   = quantile(yll,   alpha,    na.rm = TRUE),
      yll_upr   = quantile(yll, 1 - alpha,  na.rm = TRUE),
      le_m0_lwr = quantile(le_m0, alpha,    na.rm = TRUE),
      le_m0_upr = quantile(le_m0,1 - alpha, na.rm = TRUE),
      le_m1_lwr = quantile(le_m1, alpha,    na.rm = TRUE),
      le_m1_upr = quantile(le_m1,1 - alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

# 正規近似 CI を作る: bootstrap 標準誤差から信頼区間を計算する。
# Compute normal-approximation CIs: derives confidence intervals from bootstrap SEs.
yll_ci_normal <- function(point_est, boot_df, conf_level) {
  alpha <- (1 - conf_level) / 2
  z <- qnorm(1 - alpha)
  
  se_df <- boot_df |>
    group_by(age_start) |>
    summarise(
      se_yll   = sd(yll,   na.rm = TRUE),
      se_le_m0 = sd(le_m0, na.rm = TRUE),
      se_le_m1 = sd(le_m1, na.rm = TRUE),
      .groups = "drop"
    )
  
  point_est |>
    left_join(se_df, by = "age_start") |>
    transmute(
      age_start,
      se_yll, se_le_m0, se_le_m1,
      yll_lwr   = yll   - z * se_yll,
      yll_upr   = yll   + z * se_yll,
      le_m0_lwr = le_m0 - z * se_le_m0,
      le_m0_upr = le_m0 + z * se_le_m0,
      le_m1_lwr = le_m1 - z * se_le_m1,
      le_m1_upr = le_m1 + z * se_le_m1
    )
}

# 読みやすい列名に直す: 結果テーブルを用途別に分かりやすい名前へ変換する。
# Make result tables readable: renames output columns to more descriptive names.
yll_make_readable_results <- function(point_est, boot_df, summary_df, method) {
  detailed_results <- summary_df |>
    rename(
      starting_age = age_start,
      yll = yll,
      le_reference = le_m0,
      le_exposed = le_m1
    ) |>
    rename_with(~ sub("^se_yll_norm$", "yll_se_normal", .x)) |>
    rename_with(~ sub("^se_le_m0_norm$", "le_reference_se_normal", .x)) |>
    rename_with(~ sub("^se_le_m1_norm$", "le_exposed_se_normal", .x)) |>
    rename_with(~ sub("^yll_lwr_perc$", "ci_low_percentile", .x)) |>
    rename_with(~ sub("^yll_upr_perc$", "ci_high_percentile", .x)) |>
    rename_with(~ sub("^le_m0_lwr_perc$", "le_reference_ci_low_percentile", .x)) |>
    rename_with(~ sub("^le_m0_upr_perc$", "le_reference_ci_high_percentile", .x)) |>
    rename_with(~ sub("^le_m1_lwr_perc$", "le_exposed_ci_low_percentile", .x)) |>
    rename_with(~ sub("^le_m1_upr_perc$", "le_exposed_ci_high_percentile", .x)) |>
    rename_with(~ sub("^yll_lwr_norm$", "ci_low_normal", .x)) |>
    rename_with(~ sub("^yll_upr_norm$", "ci_high_normal", .x)) |>
    rename_with(~ sub("^le_m0_lwr_norm$", "le_reference_ci_low_normal", .x)) |>
    rename_with(~ sub("^le_m0_upr_norm$", "le_reference_ci_high_normal", .x)) |>
    rename_with(~ sub("^le_m1_lwr_norm$", "le_exposed_ci_low_normal", .x)) |>
    rename_with(~ sub("^le_m1_upr_norm$", "le_exposed_ci_high_normal", .x))
  
  main_results <- detailed_results |>
    select(starting_age, yll) |>
    mutate(
      ci_low = NA_real_,
      ci_high = NA_real_,
      ci_method = method
    )
  
  if (identical(method, "percentile") &&
      all(c("ci_low_percentile", "ci_high_percentile") %in% names(detailed_results))) {
    main_results <- detailed_results |>
      transmute(
        starting_age,
        yll,
        ci_low = ci_low_percentile,
        ci_high = ci_high_percentile,
        ci_method = method
      )
  } else if (identical(method, "normal") &&
             all(c("ci_low_normal", "ci_high_normal") %in% names(detailed_results))) {
    main_results <- detailed_results |>
      transmute(
        starting_age,
        yll,
        ci_low = ci_low_normal,
        ci_high = ci_high_normal,
        ci_method = method
      )
  }
  
  list(
    detailed_results = detailed_results,
    summary = main_results
  )
}

# 結果オブジェクトを作る: 主結果・詳細表・設定情報を返す。
# Build result object: returns the main table, detailed results, and metadata.
yll_build_result_object <- function(point_est, boot_df, summary_df, meta, method) {
  readable <- yll_make_readable_results(point_est, boot_df, summary_df, method = method)
  
  list(
    detailed_results = readable$detailed_results,
    summary = readable$summary,
    meta = meta
  )
}

# 既存メイン関数: ATT/ATC/ATE を指定して YLL と bootstrap CI を返す。
# Main legacy interface: estimates YLL and bootstrap CIs for ATT/ATC/ATE.
estimate_yll_gformula <- function(
    data,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    method = c("normal", "percentile"),
    show_progress = TRUE,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    estimand = c("ATT", "ATC", "ATE"),
    # 並列は呼び出し側で future::plan() を設定する前提
    use_future = TRUE
) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  
  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  
  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level
  
  set.seed(seed)
  
  ids <- unique(data[[id_var]])
  
  point_est <- etimate_yll_gformula_single(
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
    estimand = estimand
  )
  
  if (B < 1L) {
    warning("B = 0 so confidence intervals were not computed.", call. = FALSE)
    
    return(yll_build_result_object(
      point_est = point_est,
      boot_df = tibble(),
      summary_df = arrange(point_est, age_start),
      meta = list(B = B, seed = seed, conf_level = conf_level, method = method, estimand = estimand),
      method = method
    ))
  }
  
  one_boot <- function(b) {
    sampled_ids <- yll_bootstrap_ids(ids)
    boot_data <- yll_make_boot_data(data, id_var, sampled_ids)
    
    etimate_yll_gformula_single(
      data = boot_data,
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
      estimand = estimand
    ) |>
      mutate(b = b)
  }
  
  if (use_future) {
    if (show_progress) {
      handlers(global = TRUE)
      handlers("txtprogressbar")
      boot_list <- with_progress({
        p <- progressor(steps = B)
        future_lapply(seq_len(B), function(b) { p(); one_boot(b) }, future.seed = TRUE)
      })
    } else {
      boot_list <- future_lapply(seq_len(B), one_boot, future.seed = TRUE)
    }
  } else {
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)
      boot_list <- lapply(seq_len(B), function(b) { setTxtProgressBar(pb, b); one_boot(b) })
    } else {
      boot_list <- lapply(seq_len(B), one_boot)
    }
  }
  
  boot_df <- bind_rows(boot_list)
  
  summary <- point_est
  
  if (identical(method, "percentile")) {
    perc <- yll_ci_percentile(boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_perc"), -age_start)
    
    summary <- left_join(summary, perc, by = "age_start")
  }
  
  if (identical(method, "normal")) {
    norm <- yll_ci_normal(point_est, boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_norm"), -age_start)
    
    summary <- left_join(summary, norm, by = "age_start")
  }
  
  summary <- arrange(summary, age_start)
  
  yll_build_result_object(
    point_est = point_est,
    boot_df = boot_df,
    summary_df = summary,
    meta = list(B = B, seed = seed, conf_level = conf_level, method = method, estimand = estimand),
    method = method
  )
}

# 汎用介入関数: 参照介入・比較介入・対象集団を自由に指定して YLL を推定する。
# Generic intervention interface: estimates YLL under user-defined reference/exposed interventions and target population.
estimate_yll_gformula_intervention <- function(
    data,
    intervention_reference,
    intervention_exposed,
    target_population = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    method = c("normal", "percentile"),
    show_progress = TRUE,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    use_future = TRUE
) {
  method <- match.arg(method)
  
  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  
  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level
  
  set.seed(seed)
  
  ids <- unique(data[[id_var]])
  
  point_est <- estimate_yll_gformula_engine_single(
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
    intervention_reference = intervention_reference,
    intervention_exposed = intervention_exposed,
    target_population = target_population
  )
  
  if (B < 1L) {
    warning("B = 0 so confidence intervals were not computed.", call. = FALSE)
    
    return(yll_build_result_object(
      point_est = point_est,
      boot_df = tibble(),
      summary_df = arrange(point_est, age_start),
      meta = list(B = B, seed = seed, conf_level = conf_level, method = method),
      method = method
    ))
  }
  
  one_boot <- function(b) {
    sampled_ids <- yll_bootstrap_ids(ids)
    boot_data <- yll_make_boot_data(data, id_var, sampled_ids)
    
    estimate_yll_gformula_engine_single(
      data = boot_data,
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
      intervention_reference = intervention_reference,
      intervention_exposed = intervention_exposed,
      target_population = target_population
    ) |>
      mutate(b = b)
  }
  
  if (use_future) {
    if (show_progress) {
      handlers(global = TRUE)
      handlers("txtprogressbar")
      boot_list <- with_progress({
        p <- progressor(steps = B)
        future_lapply(seq_len(B), function(b) { p(); one_boot(b) }, future.seed = TRUE)
      })
    } else {
      boot_list <- future_lapply(seq_len(B), one_boot, future.seed = TRUE)
    }
  } else {
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)
      boot_list <- lapply(seq_len(B), function(b) { setTxtProgressBar(pb, b); one_boot(b) })
    } else {
      boot_list <- lapply(seq_len(B), one_boot)
    }
  }
  
  boot_df <- bind_rows(boot_list)
  
  summary <- point_est
  
  if (identical(method, "percentile")) {
    perc <- yll_ci_percentile(boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_perc"), -age_start)
    
    summary <- left_join(summary, perc, by = "age_start")
  }
  
  if (identical(method, "normal")) {
    norm <- yll_ci_normal(point_est, boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_norm"), -age_start)
    
    summary <- left_join(summary, norm, by = "age_start")
  }
  
  summary <- arrange(summary, age_start)
  
  yll_build_result_object(
    point_est = point_est,
    boot_df = boot_df,
    summary_df = summary,
    meta = list(B = B, seed = seed, conf_level = conf_level, method = method),
    method = method
  )
}

# 二値 stochastic 介入ラッパー: 2つの介入シナリオを「非曝露時」「曝露時」の曝露確率で指定して比較する。
# Binary stochastic wrapper: compares two intervention regimes defined by exposure probabilities among the currently unexposed/exposed.
estimate_yll_gformula_binary_stochastic <- function(
    data,
    prob_exposed_if_unexposed_reference,
    prob_exposed_if_exposed_reference,
    prob_exposed_if_unexposed_exposed,
    prob_exposed_if_exposed_exposed,
    target_population = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    method = c("normal", "percentile"),
    show_progress = TRUE,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    use_future = TRUE
) {
  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level
  
  intervention_reference <- yll_make_binary_stochastic_intervention(
    prob_exposed_if_unexposed = prob_exposed_if_unexposed_reference,
    prob_exposed_if_exposed = prob_exposed_if_exposed_reference,
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  
  intervention_exposed <- yll_make_binary_stochastic_intervention(
    prob_exposed_if_unexposed = prob_exposed_if_unexposed_exposed,
    prob_exposed_if_exposed = prob_exposed_if_exposed_exposed,
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  
  estimate_yll_gformula_intervention(
    data = data,
    intervention_reference = intervention_reference,
    intervention_exposed = intervention_exposed,
    target_population = target_population,
    B = B,
    seed = seed,
    conf_level = conf_level,
    method = method,
    show_progress = show_progress,
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
    use_future = use_future
  )
}

# 自然経過との比較ラッパー: binary stochastic 介入を「現状維持」と比較する。
# Natural-course comparison wrapper: compares a binary stochastic intervention against the observed baseline exposure pattern.
estimate_yll_gformula_binary_stochastic_vs_natural <- function(
    data,
    prob_exposed_if_unexposed,
    prob_exposed_if_exposed,
    target_population = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    method = c("normal", "percentile"),
    show_progress = TRUE,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    use_future = TRUE
) {
  estimate_yll_gformula_binary_stochastic(
    data = data,
    prob_exposed_if_unexposed_reference = 0,
    prob_exposed_if_exposed_reference = 1,
    prob_exposed_if_unexposed_exposed = prob_exposed_if_unexposed,
    prob_exposed_if_exposed_exposed = prob_exposed_if_exposed,
    target_population = target_population,
    B = B,
    seed = seed,
    conf_level = conf_level,
    method = method,
    show_progress = show_progress,
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
    use_future = use_future
  )
}

# ATE ラッパー: 集団全体を対象に、全員を参照曝露と比較曝露に割り付けたときの YLL を推定する。
# ATE wrapper: estimates YLL for the full population under the two exposure assignments.
estimate_yll_gformula_ate <- function(...) {
  estimate_yll_gformula(..., estimand = "ATE")
}

# ATT ラッパー: 実際に曝露群だった集団を対象に、参照曝露と比較曝露の YLL を推定する。
# ATT wrapper: estimates YLL for those actually exposed, comparing reference vs exposed assignments.
estimate_yll_gformula_att <- function(...) {
  estimate_yll_gformula(..., estimand = "ATT")
}

# ATC ラッパー: 実際に非曝露群だった集団を対象に、参照曝露と比較曝露の YLL を推定する。
# ATC wrapper: estimates YLL for those actually unexposed, comparing reference vs exposed assignments.
estimate_yll_gformula_atc <- function(...) {
  estimate_yll_gformula(..., estimand = "ATC")
}
