# Reconstruct the per-iteration conditional survival curves from the stack of
# bootstrap marginal curves.
#
# `marginal_survival_boot` is a long tibble with one row per (iteration `b`,
# age) and the two marginal arms `surv0`/`surv1`. For each iteration we slice
# out that iteration's curve, condition it on survival to `age_start`
# (i.e. divide by S_b(a_start)), and tag the result with `b` so callers can
# compute pointwise quantiles across iterations. NULL-safe so degenerate
# iterations (e.g. S_b(a_start) = 0) are silently dropped instead of poisoning
# the downstream summaries.
#' @noRd
yll_bootstrap_conditional_curves <- function(marginal_survival_boot, age_start) {
  if (is.null(marginal_survival_boot) || nrow(marginal_survival_boot) == 0) {
    return(NULL)
  }

  bs <- unique(marginal_survival_boot[["b"]])
  out_list <- lapply(bs, function(bi) {
    one_b <- marginal_survival_boot[marginal_survival_boot[["b"]] == bi, , drop = FALSE]
    cs <- yll_conditional_survival_from_age(one_b, age_start)
    if (is.null(cs)) return(NULL)
    cs[["b"]] <- bi
    cs
  })
  out_list <- out_list[!vapply(out_list, is.null, logical(1))]
  if (length(out_list) == 0) return(NULL)
  bind_rows(out_list)
}

# Build the conditional-survival table the plotting helper consumes.
#
# Combines:
#   * the point estimate `S(y | a_start)` for each arm (computed from the
#     point-estimate marginal curves), and
#   * pointwise confidence bounds at every age, taken as the lower/upper
#     `(1 - conf_level)/2` quantiles of the bootstrap conditional curves.
#
# When no bootstrap curves are available (e.g. `B = 0`) we return just the
# point estimate columns, and the calling plot function will skip the ribbon.
#' @noRd
yll_conditional_curves_with_ci <- function(marginal_survival_point,
                                           marginal_survival_boot,
                                           age_start,
                                           conf_level = 0.95) {
  if (is.null(marginal_survival_point)) {
    stop("Result has no marginal_survival_point. Re-run estimation to enable plotting.", call. = FALSE)
  }

  point_cond <- yll_conditional_survival_from_age(marginal_survival_point, age_start)
  if (is.null(point_cond)) {
    stop("age_start = ", age_start, " is outside the available range or marginal survival is 0.", call. = FALSE)
  }

  out <- tibble(
    age_temp     = point_cond$age_temp,
    surv_cond_g0 = point_cond$surv_cond_g0,
    surv_cond_g1 = point_cond$surv_cond_g1
  )

  boot_cond <- yll_bootstrap_conditional_curves(marginal_survival_boot, age_start)
  if (!is.null(boot_cond) && nrow(boot_cond) > 0) {
    alpha <- (1 - conf_level) / 2
    ci_df <- boot_cond |>
      group_by(age_temp) |>
      summarise(
        surv_cond_g0_low  = quantile(surv_cond_g0, alpha,    na.rm = TRUE),
        surv_cond_g0_high = quantile(surv_cond_g0, 1 - alpha, na.rm = TRUE),
        surv_cond_g1_low  = quantile(surv_cond_g1, alpha,    na.rm = TRUE),
        surv_cond_g1_high = quantile(surv_cond_g1, 1 - alpha, na.rm = TRUE),
        .groups = "drop"
      )
    out <- left_join(out, ci_df, by = "age_temp")
  }

  out
}

# Build the marginal-survival table for `plot_marginal_survival()`.
#
# Mirrors `yll_conditional_curves_with_ci()` but on the *unconditional*
# population-marginal curves S(t) — i.e. before conditioning on survival to
# any chosen starting age. Pointwise CI bounds at each age are again the
# lower/upper `(1 - conf_level)/2` quantiles across bootstrap iterations.
#' @noRd
yll_marginal_curves_with_ci <- function(marginal_survival_point,
                                        marginal_survival_boot,
                                        conf_level = 0.95) {
  if (is.null(marginal_survival_point)) {
    stop("Result has no marginal_survival_point. Re-run estimation to enable plotting.", call. = FALSE)
  }

  out <- tibble(
    age_temp = marginal_survival_point$age_temp,
    surv0    = marginal_survival_point$surv0,
    surv1    = marginal_survival_point$surv1
  )

  if (!is.null(marginal_survival_boot) && nrow(marginal_survival_boot) > 0) {
    alpha <- (1 - conf_level) / 2
    ci_df <- marginal_survival_boot |>
      group_by(age_temp) |>
      summarise(
        surv0_low  = quantile(surv0, alpha,    na.rm = TRUE),
        surv0_high = quantile(surv0, 1 - alpha, na.rm = TRUE),
        surv1_low  = quantile(surv1, alpha,    na.rm = TRUE),
        surv1_high = quantile(surv1, 1 - alpha, na.rm = TRUE),
        .groups = "drop"
      )
    out <- left_join(out, ci_df, by = "age_temp")
  }

  out
}

#' Plot the conditional survival probability from a starting age
#'
#' Draws the conditional survival curves \eqn{S(y \mid a_{\text{start}})} for the
#' two intervention arms, optionally with a pointwise bootstrap confidence
#' band.
#'
#' @param res A result object returned by [estimate_yll_gformula()] (or one of
#'   its variants).
#' @param age_start Numeric. The starting age \eqn{a_{\text{start}}} used to
#'   condition the survival curve.
#' @param conf_band Logical. If `TRUE` (default) and bootstrap curves are
#'   available, draws a pointwise confidence band around each curve.
#' @param reference_label,exposed_label Character labels for the two
#'   intervention arms shown in the legend.
#' @param conf_level Confidence level for the pointwise band. Defaults to the
#'   value stored in `res$meta$conf_level`, falling back to `0.95`.
#'
#' @return A `ggplot` object.
#' @export
plot_conditional_survival <- function(res,
                                      age_start,
                                      conf_band = TRUE,
                                      reference_label = "Reference",
                                      exposed_label   = "Exposed",
                                      conf_level      = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_conditional_survival().", call. = FALSE)
  }

  if (is.null(conf_level)) {
    conf_level <- if (!is.null(res$meta$conf_level)) res$meta$conf_level else 0.95
  }

  cs <- yll_conditional_curves_with_ci(
    marginal_survival_point = res$marginal_survival_point,
    marginal_survival_boot  = res$marginal_survival_boot,
    age_start               = age_start,
    conf_level              = conf_level
  )

  has_ci <- all(c("surv_cond_g0_low", "surv_cond_g0_high",
                  "surv_cond_g1_low", "surv_cond_g1_high") %in% names(cs))

  long <- bind_rows(
    tibble(
      age_temp = cs$age_temp,
      surv     = cs$surv_cond_g0,
      ci_low   = if (has_ci) cs$surv_cond_g0_low  else NA_real_,
      ci_high  = if (has_ci) cs$surv_cond_g0_high else NA_real_,
      arm      = reference_label
    ),
    tibble(
      age_temp = cs$age_temp,
      surv     = cs$surv_cond_g1,
      ci_low   = if (has_ci) cs$surv_cond_g1_low  else NA_real_,
      ci_high  = if (has_ci) cs$surv_cond_g1_high else NA_real_,
      arm      = exposed_label
    )
  )
  long$arm <- factor(long$arm, levels = c(reference_label, exposed_label))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = age_temp, y = surv,
                                          colour = arm, fill = arm)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x      = "Age (years)",
      y      = sprintf("Conditional survival probability (from age %g)", age_start),
      colour = "Intervention",
      fill   = "Intervention"
    ) +
    ggplot2::theme_minimal()

  if (conf_band && has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}

#' Plot YLL across starting ages
#'
#' Draws estimated years of life lost (YLL) as a function of the starting age,
#' with an optional confidence band when bootstrap CIs are available in the
#' result object.
#'
#' @param res A result object returned by [estimate_yll_gformula()] (or one of
#'   its variants).
#' @param conf_band Logical. If `TRUE` (default) and `ci_low`/`ci_high` are
#'   present in `res$summary`, draws a confidence band.
#'
#' @return A `ggplot` object.
#' @export
plot_yll <- function(res, conf_band = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_yll().", call. = FALSE)
  }

  s <- res$summary

  p <- ggplot2::ggplot(s, ggplot2::aes(x = starting_age, y = yll)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Starting age (years)",
      y = "Years of life lost (years)"
    ) +
    ggplot2::theme_minimal()

  if (conf_band && all(c("ci_low", "ci_high") %in% names(s)) &&
      any(is.finite(s$ci_low))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}

#' Plot marginal survival curves for the two intervention arms
#'
#' Draws the population-level marginal survival curves \eqn{S(t)} for the two
#' intervention arms, optionally with pointwise bootstrap confidence bands.
#'
#' @inheritParams plot_conditional_survival
#'
#' @return A `ggplot` object.
#' @export
plot_marginal_survival <- function(res,
                                   conf_band = TRUE,
                                   reference_label = "Reference",
                                   exposed_label   = "Exposed",
                                   conf_level      = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_marginal_survival().", call. = FALSE)
  }

  if (is.null(conf_level)) {
    conf_level <- if (!is.null(res$meta$conf_level)) res$meta$conf_level else 0.95
  }

  ms <- yll_marginal_curves_with_ci(
    marginal_survival_point = res$marginal_survival_point,
    marginal_survival_boot  = res$marginal_survival_boot,
    conf_level              = conf_level
  )

  has_ci <- all(c("surv0_low", "surv0_high", "surv1_low", "surv1_high") %in% names(ms))

  long <- bind_rows(
    tibble(
      age_temp = ms$age_temp,
      surv     = ms$surv0,
      ci_low   = if (has_ci) ms$surv0_low  else NA_real_,
      ci_high  = if (has_ci) ms$surv0_high else NA_real_,
      arm      = reference_label
    ),
    tibble(
      age_temp = ms$age_temp,
      surv     = ms$surv1,
      ci_low   = if (has_ci) ms$surv1_low  else NA_real_,
      ci_high  = if (has_ci) ms$surv1_high else NA_real_,
      arm      = exposed_label
    )
  )
  long$arm <- factor(long$arm, levels = c(reference_label, exposed_label))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = age_temp, y = surv,
                                          colour = arm, fill = arm)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x      = "Age (years)",
      y      = "Marginal survival probability",
      colour = "Intervention",
      fill   = "Intervention"
    ) +
    ggplot2::theme_minimal()

  if (conf_band && has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}
