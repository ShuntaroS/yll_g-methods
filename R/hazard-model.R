# カット点を作る: 離散時間化のための分割点を追跡時間から作成する。
# Build cut points: creates split points for discretizing follow-up time.
#' @noRd
yll_make_cut_points <- function(data, time_var, by = 1L) {
  max_t <- ceiling(max(data[[time_var]], na.rm = TRUE))
  seq.int(from = 0L, to = max_t, by = by)
}

# 生存データをロング形式へ変換する: `survSplit()` で時間区間ごとのデータに展開する。
# Split survival data: expands person-level data into interval-level rows with `survSplit()`.
#' @noRd
yll_survsplit <- function(data, time_var, event_var, cut_points) {
  surv_formula <- stats::as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ."))
  survival::survSplit(
    formula  = surv_formula,
    data     = data,
    cut      = cut_points[-1],
    episode  = "tgroup"
  )
}

# ノット候補を作る: 年齢変数の分位点から spline のノットを作成する。
# Build spline knots: uses quantiles of age to define spline knots.
#' @noRd
yll_knot_quantiles <- function(x, probs = c(0, 0.25, 0.50, 0.75, 1)) {
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE))
}

# 交絡因子の右辺項を作る: モデル式に入れる共変量部分を文字列で返す。
# Build confounder RHS: returns the covariate part of the model formula.
#' @noRd
yll_rhs_confounders <- function(confounders_baseline) {
  if (is.null(confounders_baseline) || length(confounders_baseline) == 0) {
    "1"
  } else {
    paste(confounders_baseline, collapse = " + ")
  }
}

# 離散時間ハザードモデルを当てる: 年齢 spline と曝露の交互作用を含む GLM を推定する。
# Fit discrete-time hazard model: estimates a GLM with age splines and exposure interaction.
#' @noRd
yll_fit_discrete_hazard_glm <- function(
    df_long,
    event_var,
    confounders_baseline = NULL,
    knot_p
) {
  rhs <- yll_rhs_confounders(confounders_baseline)

  form <- stats::as.formula(
    paste0(
      event_var,
      " ~ ",
      "expo",
      " + Epi::Ns(age_temp, knots = knot_p) + ",
      "expo * Epi::Ns(age_temp, knots = knot_p) + ",
      rhs
    )
  )

  stats::glm(formula = form, family = stats::binomial(), data = df_long)
}

# ハザードを予測する: 学習済みモデルから各行の死亡確率を返す。
# Predict hazards: returns row-wise event probabilities from the fitted model.
#' @noRd
yll_predict_hazard <- function(fit, newdata) {
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}
