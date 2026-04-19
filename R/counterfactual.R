# 反実仮想データを展開する: 各個体について `age_temp_start..age_temp_end` の予測用行を作る。
# Expand counterfactual data: creates one prediction row per individual for each age in `age_temp_start..age_temp_end`.
#' @noRd
yll_expand_counterfactual_data <- function(
    data,
    id_var,
    exposure_var,
    age_temp_start,
    age_temp_end
) {
  if (age_temp_end < age_temp_start) {
    stop("`age_temp_end` must be >= `age_temp_start`.", call. = FALSE)
  }
  n_ages <- as.integer(age_temp_end - age_temp_start + 1L)

  out <- data |>
    tidyr::uncount(weights = n_ages, .remove = FALSE) |>
    dplyr::group_by(.data[[id_var]]) |>
    dplyr::mutate(age_temp = age_temp_start + dplyr::row_number() - 1L) |>
    dplyr::ungroup() |>
    dplyr::rename(expo_original = dplyr::all_of(exposure_var))

  out[["expo"]] <- out$expo_original
  out
}

# 個体別生存を作る: 参照曝露と曝露群の下での生存曲線を個体ごとに計算する。
# Build individual survival curves: computes per-subject survival under the reference and exposed levels.
#' @noRd
yll_compute_individual_survival_curves <- function(df, id_var,
                                                   hazard0_var = "hazard0",
                                                   hazard1_var = "hazard1") {
  df |>
    arrange(.data[[id_var]], age_temp) |>
    group_by(.data[[id_var]]) |>
    mutate(
      .surv0_after_interval = cumprod(1 - .data[[hazard0_var]]),
      .surv1_after_interval = cumprod(1 - .data[[hazard1_var]]),
      surv0 = lag(.surv0_after_interval, default = 1),
      surv1 = lag(.surv1_after_interval, default = 1)
    ) |>
    ungroup()
}

# 介入下の平均生存を作る: 個体別の曝露確率で2本の生存曲線を混合して平均する。
# Build intervention survival: mixes the two survival curves using subject-specific exposure probabilities.
#' @noRd
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
#' @noRd
yll_mean_survival_by_age <- function(df, id_var, hazard_var, surv_name = "surv") {
  df |>
    arrange(.data[[id_var]], age_temp) |>
    group_by(.data[[id_var]]) |>
    mutate(
      .surv_after_interval = cumprod(1 - .data[[hazard_var]]),
      .surv_at_age = lag(.surv_after_interval, default = 1)
    ) |>
    ungroup() |>
    group_by(age_temp) |>
    summarise(!!surv_name := mean(.surv_at_age), .groups = "drop")
}

# 条件付き生存を作る: 開始年齢で標準化した生存曲線に変換する。
# Build conditional survival: rescales survival curves to start at a chosen age.
#' @noRd
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
# `integration` selects the discretisation rule:
#   - "left_rectangle": sum_i (a_{i+1} - a_i) * S(a_i). 各区間の開始時点の生存をその区間の代表値とする。
#     離散時間 hazard モデルの自然な対応で、デフォルトはこちら。
#   - "trapezoidal":    sum_i (a_{i+1} - a_i) * (S(a_i) + S(a_{i+1})) / 2.
#     生存曲線を区間内で線形補間した近似で、滑らかな曲線では精度がやや高い。
#' @noRd
yll_integrate_le <- function(age, surv,
                             integration = c("left_rectangle", "trapezoidal")) {
  integration <- match.arg(integration)

  if (length(age) <= 1) {
    return(0)
  }

  ord <- order(age)
  age <- age[ord]
  surv <- surv[ord]

  d_age <- diff(age)

  if (identical(integration, "left_rectangle")) {
    return(sum(d_age * head(surv, -1)))
  }

  sum(d_age * (head(surv, -1) + tail(surv, -1)) / 2)
}
