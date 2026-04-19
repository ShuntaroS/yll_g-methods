# Translate one of the three named estimands (ATE / ATT / ATC) into a
# population-selection rule used at averaging time.
#
#   * ATE: average over the entire counterfactual population => no rule.
#   * ATT: average only over subjects whose *observed* exposure was the
#          exposed level (treatment effect on the treated).
#   * ATC: average only over subjects whose *observed* exposure was the
#          reference level (treatment effect on the controls).
#
# Note that this rule is applied *after* counterfactual survival curves have
# been built for everyone — never to subset the model-fitting data. That
# matters: the hazard model is always fit on the full sample, and the
# estimand only changes which subjects' counterfactual contrasts we average.
#' @noRd
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

# Resolve a target-population specification into the logical row-mask that
# the engine will use to subset the counterfactual data before averaging.
#
# `target_population` may be either a logical vector (already pre-computed)
# or a function `f(data)` returning a logical vector. We validate length and
# the absence of NAs because silently dropping rows here would change the
# estimand without warning.
#' @noRd
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

# Convert a user-supplied `intervention_*` argument into a numeric vector of
# subject-specific exposure probabilities `P(A* = exposed)`. The accepted
# input forms are:
#
#   * numeric scalar in [0, 1] — same probability for everyone (broadcast).
#   * numeric vector of length nrow(data) — already subject-specific.
#   * logical scalar/vector — coerced via `as.numeric`.
#   * a value matching one of the two exposure levels — interpreted as a
#     deterministic intervention assigning that level to everyone.
#   * a function `f(data)` returning any of the above.
#
# The returned vector is what the engine multiplies the per-arm survival
# curves by to form the marginal mixture S^* = (1 - p) * S_ref + p * S_exp.
#' @noRd
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

#' Build a binary stochastic intervention
#'
#' Creates an intervention function suitable for the
#' `intervention_reference` / `intervention_exposed` arguments of
#' [estimate_yll_gformula_intervention()]. For each subject the returned
#' function reports `P(A* = exposed_level)` as a function of the subject's
#' originally observed exposure status.
#'
#' @param prob_exposed_if_unexposed Probability of being assigned the exposed
#'   level among subjects whose observed exposure was the reference level.
#' @param prob_exposed_if_exposed Probability of being assigned the exposed
#'   level among subjects whose observed exposure was the exposed level.
#' @param reference_level,exposed_level Values labelling the two exposure levels
#'   in the original data.
#'
#' @return A function `f(data)` returning a numeric vector of length
#'   `nrow(data)` whose entries are in `[0, 1]`.
#'
#' @export
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
