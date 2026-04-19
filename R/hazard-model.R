# Build the cut points used to discretise follow-up time into one-year (or
# `by`-year) intervals. We always start the grid at 0 and extend it to
# ceiling(max(time)) so that every observed event/censoring time is covered.
#' @noRd
yll_make_cut_points <- function(data, time_var, by = 1L) {
  max_t <- ceiling(max(data[[time_var]], na.rm = TRUE))
  seq.int(from = 0L, to = max_t, by = by)
}

# Convert person-level (id, time, event) data into person-period (long) data,
# one row per individual per follow-up interval. This is the standard input
# shape for a discrete-time pooled logistic hazard model: each row carries the
# event indicator for that specific interval, and we'll later regress it on
# the time-varying age via splines.
#
# `survival::survSplit()` does the heavy lifting; we drop the very first cut
# point (0) because it is the start of follow-up, not a split.
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

# Choose knot positions for the natural cubic spline on attained age. Using
# data-driven quantiles (rather than a fixed grid) keeps the spline well-
# supported even when the age distribution is skewed. Defaults are min, the
# three quartiles, and max — the conventional 5-knot setup for `Epi::Ns`.
#' @noRd
yll_knot_quantiles <- function(x, probs = c(0, 0.25, 0.50, 0.75, 1)) {
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE))
}

# Build the right-hand side string for the baseline-confounder portion of the
# hazard formula. Returning "1" for the empty case keeps the formula valid
# (`event ~ expo + Ns(...) + expo * Ns(...) + 1`) and means the rest of the
# pipeline doesn't need to special-case "no confounders".
#' @noRd
yll_rhs_confounders <- function(confounders_baseline) {
  if (is.null(confounders_baseline) || length(confounders_baseline) == 0) {
    "1"
  } else {
    paste(confounders_baseline, collapse = " + ")
  }
}

# Fit the discrete-time hazard model used throughout the package: a pooled
# logistic regression of the per-interval event indicator on
#
#   * the binary exposure `expo`,
#   * a natural cubic spline in attained age (`age_temp`),
#   * the spline-by-exposure interaction (so the hazard ratio is allowed to
#     vary with age), and
#   * any user-supplied baseline confounders.
#
# The interaction term is what makes the counterfactual contrast non-trivial:
# without it the exposure effect would collapse to a single hazard ratio.
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

# Predict per-row event probabilities from the fitted hazard model. We force
# `type = "response"` so the result is on the probability scale (these will
# later be combined into survival curves via cumulative-product).
#' @noRd
yll_predict_hazard <- function(fit, newdata) {
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}
