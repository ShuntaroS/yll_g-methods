# Check that every variable referenced by the user-facing API is actually
# present in `data`. Failing fast here gives a much clearer error message than
# the cryptic `object 'foo' not found` that would otherwise come from inside
# the modelling pipeline.
#' @noRd
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

# Warn the user if they passed the entry-age column inside
# `confounders_baseline`. The hazard model already includes a flexible spline
# in attained age (= entry age + time-on-study), so adding entry age as a
# baseline confounder usually causes two problems:
#
#   (i) it splits the age effect into a redundant entry-age term and a hidden
#       time-on-study term, which makes the contrast harder to interpret; and
#  (ii) it forces the hazard model to extrapolate when we predict
#       counterfactual hazards at age_temp < entry age (the counterfactual
#       expansion creates rows for every age in [age_start, age_end] for every
#       subject, including ages before they actually entered the study).
#
# We only warn rather than error, because the user may genuinely want a
# cohort / time-on-study effect.
#' @noRd
yll_warn_entry_age_in_confounders <- function(age_at_entry_var, confounders_baseline) {
  if (is.null(confounders_baseline) || length(confounders_baseline) == 0) {
    return(invisible(NULL))
  }
  if (age_at_entry_var %in% confounders_baseline) {
    warning(
      "`", age_at_entry_var, "` is included in `confounders_baseline`. ",
      "Because attained age (= `", age_at_entry_var, "` + time-on-study) is already ",
      "modelled via splines, including entry age typically (i) double-models the age ",
      "effect via a hidden time-on-study term, and (ii) forces extrapolation in ",
      "counterfactual hazard prediction at age_temp < `", age_at_entry_var, "`. ",
      "Consider removing it unless you intentionally want to model time-on-study or ",
      "cohort effects.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Decide which observed value is the "reference" level (= unexposed) and which
# is the "exposed" level for a binary exposure column.
#
# Behaviour:
#   * If the user supplies both `reference_level` and `exposed_level` we just
#     validate them.
#   * If either is NULL we infer them: for a factor we use the natural
#     `levels()` order; for a non-factor we use the order of unique observed
#     values. The first becomes the reference, the second the exposed.
#
# The package currently only supports binary exposures, so anything other than
# exactly two distinct non-NA values is an error.
#' @noRd
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
