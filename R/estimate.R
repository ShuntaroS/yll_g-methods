#' Estimate YLL via the g-formula (ATE / ATT / ATC)
#'
#' Fits a discrete-time hazard model on the attained-age time scale, builds
#' counterfactual survival curves under "all reference" vs. "all exposed"
#' interventions, and reports Years of Life Lost (YLL) and counterfactual life
#' expectancy (LE) at user-chosen starting ages, with bootstrap confidence
#' intervals.
#'
#' This function is the main legacy interface. For arbitrary interventions or
#' target populations, use [estimate_yll_gformula_intervention()].
#'
#' @param data A `data.frame` with one row per individual, containing the
#'   variables specified by `id_var`, `time_var`, `event_var`, `exposure_var`,
#'   `age_at_entry_var`, and any `confounders_baseline`.
#' @param B Number of bootstrap iterations. Set to `0` to skip bootstrapping
#'   (returns point estimates only).
#' @param seed Integer seed for reproducibility.
#' @param conf_level Confidence level for the intervals (e.g. `0.95`).
#' @param method Either `"normal"` (Wald-type interval using bootstrap SE) or
#'   `"percentile"` (bootstrap percentile interval).
#' @param show_progress Show a progress bar during bootstrap.
#' @param id_var,time_var,event_var,exposure_var Column names in `data`.
#' @param reference_level,exposed_level Values labelling the two exposure
#'   levels. If `NULL` they are inferred from the data.
#' @param age_at_entry_var Column name for age at study entry.
#' @param age_start,age_end,age_interval Starting ages for which YLL is
#'   reported, given as `seq(age_start, age_end, age_interval)`.
#' @param confounders_baseline Character vector of baseline confounder column
#'   names. Including the entry-age column here is discouraged and triggers a
#'   warning, since attained age is already modelled via splines.
#' @param estimand One of `"ATT"`, `"ATC"`, or `"ATE"`.
#' @param integration Discretisation rule for the area under the conditional
#'   survival curve: `"left_rectangle"` (default; matches the discrete-time
#'   hazard model) or `"trapezoidal"`.
#' @param use_future If `TRUE`, parallelise the bootstrap with
#'   [future.apply::future_lapply()]. The caller is expected to set a parallel
#'   plan (e.g. `future::plan(future::multisession())`).
#'
#' @return A list with elements
#'   \describe{
#'     \item{`summary`}{Point estimates and CI for YLL at each starting age.}
#'     \item{`detailed_results`}{Same as `summary` plus LE under each arm and
#'       all CI / SE columns.}
#'     \item{`meta`}{Run settings (`B`, `seed`, `conf_level`, `method`,
#'       `estimand`, `age_start`, `age_end`, `age_interval`, `integration`).}
#'     \item{`marginal_survival_point`}{Population-marginal counterfactual
#'       survival curves (point estimate).}
#'     \item{`marginal_survival_boot`}{Same curves across bootstrap
#'       iterations (used by the plotting helpers).}
#'   }
#'
#' @seealso [estimate_yll_gformula_intervention()],
#'   [estimate_yll_gformula_binary_stochastic_vs_natural()], [plot_yll()].
#'
#' @examples
#' \donttest{
#' data(yll_toy, package = "YLLgmethods")
#' res <- estimate_yll_gformula_ate(
#'   data                 = yll_toy,
#'   B                    = 0,
#'   show_progress        = FALSE,
#'   time_var             = "period",
#'   event_var            = "event",
#'   exposure_var         = "hypertension",
#'   reference_level      = "No",
#'   exposed_level        = "Yes",
#'   age_at_entry_var     = "age",
#'   age_start            = 50,
#'   age_end              = 80,
#'   age_interval         = 10,
#'   confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
#'   use_future           = FALSE
#' )
#' res$summary
#' }
#'
#' @export
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
    integration = c("left_rectangle", "trapezoidal"),
    use_future = TRUE
) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  integration <- match.arg(integration)

  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  yll_warn_entry_age_in_confounders(age_at_entry_var, confounders_baseline)

  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level

  set.seed(seed)

  ids <- unique(data[[id_var]])

  point_est <- estimate_yll_gformula_single(
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
    estimand = estimand,
    integration = integration
  )
  marginal_survival_point <- attr(point_est, "marginal_curves")

  if (B < 1L) {
    warning("B = 0 so confidence intervals were not computed.", call. = FALSE)

    return(yll_build_result_object(
      point_est = point_est,
      boot_df = tibble(),
      summary_df = arrange(point_est, age_start),
      meta = list(
        B = B, seed = seed, conf_level = conf_level, method = method,
        estimand = estimand, age_start = age_start, age_end = age_end,
        age_interval = age_interval, integration = integration
      ),
      method = method,
      marginal_survival_point = marginal_survival_point,
      marginal_survival_boot  = NULL
    ))
  }

  one_boot <- function(b) {
    sampled_ids <- yll_bootstrap_ids(ids)
    boot_data <- yll_make_boot_data(data, id_var, sampled_ids)

    res_b <- estimate_yll_gformula_single(
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
      estimand = estimand,
      integration = integration
    )
    yll_one_boot_result(res_b, b = b)
  }

  boot_results <- yll_run_bootstrap(B, one_boot, use_future, show_progress)

  boot_df <- bind_rows(lapply(boot_results, `[[`, "yll"))
  marginal_survival_boot <- yll_collect_bootstrap_curves(boot_results)

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
    meta = list(
      B = B, seed = seed, conf_level = conf_level, method = method,
      estimand = estimand, age_start = age_start, age_end = age_end,
      age_interval = age_interval, integration = integration
    ),
    method = method,
    marginal_survival_point = marginal_survival_point,
    marginal_survival_boot  = marginal_survival_boot
  )
}

# bootstrap 反復を実行する: 並列・プログレスバーを切り替えて one_boot を B 回回す。
# Run bootstrap iterations: dispatches one_boot B times with optional parallelism / progress.
#' @noRd
yll_run_bootstrap <- function(B, one_boot, use_future, show_progress) {
  if (use_future) {
    if (show_progress) {
      handlers(global = TRUE)
      handlers("txtprogressbar")
      return(with_progress({
        p <- progressor(steps = B)
        future_lapply(seq_len(B), function(b) { p(); one_boot(b) }, future.seed = TRUE)
      }))
    }
    return(future_lapply(seq_len(B), one_boot, future.seed = TRUE))
  }

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    on.exit(close(pb), add = TRUE)
    return(lapply(seq_len(B), function(b) { setTxtProgressBar(pb, b); one_boot(b) }))
  }

  lapply(seq_len(B), one_boot)
}

#' Estimate YLL under arbitrary interventions and target population
#'
#' Generic g-formula interface. The two intervention scenarios and the target
#' population are user-defined; the function returns the YLL between them with
#' bootstrap CIs.
#'
#' @inheritParams estimate_yll_gformula
#' @param intervention_reference,intervention_exposed Each is one of:
#'   * a single exposure value (e.g. `"No"`),
#'   * a numeric scalar in `[0, 1]` interpreted as the marginal probability of
#'     being assigned the exposed level,
#'   * a numeric / logical vector of length `nrow(data)` of subject-specific
#'     probabilities, or
#'   * a function `f(data)` returning either of the above (see
#'     [yll_make_binary_stochastic_intervention()]).
#' @param target_population Either `NULL` (use the full population) or a
#'   function `f(data)` that returns a logical vector marking the subjects to
#'   average over. Subjects can be selected on `data$expo_original`.
#'
#' @return The same list structure as [estimate_yll_gformula()].
#'
#' @seealso [yll_make_binary_stochastic_intervention()],
#'   [estimate_yll_gformula_binary_stochastic_vs_natural()].
#'
#' @export
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
    integration = c("left_rectangle", "trapezoidal"),
    use_future = TRUE
) {
  method <- match.arg(method)
  integration <- match.arg(integration)

  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  yll_warn_entry_age_in_confounders(age_at_entry_var, confounders_baseline)

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
    target_population = target_population,
    integration = integration
  )
  marginal_survival_point <- attr(point_est, "marginal_curves")

  if (B < 1L) {
    warning("B = 0 so confidence intervals were not computed.", call. = FALSE)

    return(yll_build_result_object(
      point_est = point_est,
      boot_df = tibble(),
      summary_df = arrange(point_est, age_start),
      meta = list(
        B = B, seed = seed, conf_level = conf_level, method = method,
        age_start = age_start, age_end = age_end,
        age_interval = age_interval, integration = integration
      ),
      method = method,
      marginal_survival_point = marginal_survival_point,
      marginal_survival_boot  = NULL
    ))
  }

  one_boot <- function(b) {
    sampled_ids <- yll_bootstrap_ids(ids)
    boot_data <- yll_make_boot_data(data, id_var, sampled_ids)

    res_b <- estimate_yll_gformula_engine_single(
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
      target_population = target_population,
      integration = integration
    )
    yll_one_boot_result(res_b, b = b)
  }

  boot_results <- yll_run_bootstrap(B, one_boot, use_future, show_progress)

  boot_df <- bind_rows(lapply(boot_results, `[[`, "yll"))
  marginal_survival_boot <- yll_collect_bootstrap_curves(boot_results)

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
    meta = list(
      B = B, seed = seed, conf_level = conf_level, method = method,
      age_start = age_start, age_end = age_end,
      age_interval = age_interval, integration = integration
    ),
    method = method,
    marginal_survival_point = marginal_survival_point,
    marginal_survival_boot  = marginal_survival_boot
  )
}

#' Compare two binary stochastic interventions
#'
#' Builds the two intervention specifications from the four conditional
#' probabilities `P(A* = exposed | A_obs)` and forwards them to
#' [estimate_yll_gformula_intervention()].
#'
#' @inheritParams estimate_yll_gformula_intervention
#' @param prob_exposed_if_unexposed_reference,prob_exposed_if_exposed_reference
#'   Probabilities defining the *reference* intervention.
#' @param prob_exposed_if_unexposed_exposed,prob_exposed_if_exposed_exposed
#'   Probabilities defining the *exposed* intervention.
#'
#' @return Same list structure as [estimate_yll_gformula()].
#' @export
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
    integration = c("left_rectangle", "trapezoidal"),
    use_future = TRUE
) {
  integration <- match.arg(integration)

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
    integration = integration,
    use_future = use_future
  )
}

#' Compare a binary stochastic intervention against the natural course
#'
#' Convenience wrapper around [estimate_yll_gformula_binary_stochastic()] in
#' which the *reference* scenario is fixed to the natural course (every subject
#' keeps their observed exposure).
#'
#' @inheritParams estimate_yll_gformula_binary_stochastic
#' @param prob_exposed_if_unexposed,prob_exposed_if_exposed Probabilities
#'   defining the comparison intervention.
#'
#' @return Same list structure as [estimate_yll_gformula()].
#' @export
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
    integration = c("left_rectangle", "trapezoidal"),
    use_future = TRUE
) {
  integration <- match.arg(integration)

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
    integration = integration,
    use_future = use_future
  )
}

#' Estimand-specific wrappers
#'
#' Thin wrappers around [estimate_yll_gformula()] that pre-fill the
#' `estimand` argument.
#'
#' @param ... Arguments forwarded to [estimate_yll_gformula()].
#' @return Same list structure as [estimate_yll_gformula()].
#' @name estimate_yll_gformula_estimand_wrappers
NULL

#' @rdname estimate_yll_gformula_estimand_wrappers
#' @export
estimate_yll_gformula_ate <- function(...) {
  estimate_yll_gformula(..., estimand = "ATE")
}

#' @rdname estimate_yll_gformula_estimand_wrappers
#' @export
estimate_yll_gformula_att <- function(...) {
  estimate_yll_gformula(..., estimand = "ATT")
}

#' @rdname estimate_yll_gformula_estimand_wrappers
#' @export
estimate_yll_gformula_atc <- function(...) {
  estimate_yll_gformula(..., estimand = "ATC")
}
