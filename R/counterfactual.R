# Build the prediction grid used by the counterfactual step. For every
# subject we generate one row per integer age in [age_temp_start,
# age_temp_end]; the hazard model will then be evaluated at each of these
# rows under both intervention arms.
#
# We deliberately *do not* expand the full lifetime: shrinking the grid to
# [age_start, age_end] is mathematically harmless because the conditional
# survival S(y | a_start) cancels every hazard at ages below a_start, and it
# (a) cuts compute proportionally and (b) avoids extrapolating the hazard
# model into ages outside the support of the data.
#
# `expo_original` preserves the *observed* exposure status (used downstream
# by ATT/ATC selection rules and stochastic interventions); `expo` is set to
# the original value here and will be overwritten by the engine when it
# evaluates the two arms.
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

# Convert per-interval hazards into per-subject survival curves under the two
# intervention arms.
#
# Discrete-time identity:
#   S(t) = prod_{u < t} (1 - h(u))
#
# In code: `cumprod(1 - hazard)` gives S(t+1) at row t (survival *after* the
# interval ends). To recover S(t) — survival *at the start* of the interval,
# which is what the integration step expects — we lag by one and seed the
# first interval at S = 1.
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

# Aggregate the two per-subject survival curves into a single population-
# marginal survival curve under a *stochastic* intervention.
#
# For each subject the intervention specifies P(A* = exposed) = `p_exposed`,
# so their counterfactual survival mixes the two arms:
#
#   S_i^*(t) = (1 - p_i) * S_i^{A=ref}(t) + p_i * S_i^{A=exp}(t)
#
# We then average the per-subject mixture across subjects within each age to
# get the marginal curve. Determinstic interventions are just the special
# case `p_i in {0, 1}`.
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

# Same lag-by-one cumulative-product as `yll_compute_individual_survival_curves`
# but for a *single* hazard column, then averaged across subjects at each age.
# Used by helpers that need a single observed-distribution survival curve
# rather than the two intervention arms.
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

# Convert a marginal survival curve S(t) into the *conditional* survival
# curve given survival to a chosen starting age:
#
#   S(t | a_start) = S(t) / S(a_start) for t >= a_start.
#
# Returning NULL when S(a_start) is zero (or absent) lets the caller skip
# starting ages where the curve is undefined, instead of producing NaNs.
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

# Numerically integrate a conditional survival curve to get residual life
# expectancy (LE = area under the curve).
#
# Two discretisation rules are supported:
#
#   * "left_rectangle" (default):
#       LE ≈ sum_i (a_{i+1} - a_i) * S(a_i)
#     i.e. the survival at the *start* of each interval represents the whole
#     interval. This is the natural companion to the discrete-time hazard
#     model, where S(a) is interpreted as survival up to the start of
#     interval [a, a+1).
#
#   * "trapezoidal":
#       LE ≈ sum_i (a_{i+1} - a_i) * (S(a_i) + S(a_{i+1})) / 2
#     i.e. linearly interpolate the survival curve within each interval.
#     This is slightly more accurate when the curve is smooth but mixes a
#     continuous-time intuition into a discrete-time model.
#
# A length-<=1 `age` vector is treated as zero area, which conveniently makes
# the function safe to call on degenerate cases without special-casing them
# at the call site.
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
