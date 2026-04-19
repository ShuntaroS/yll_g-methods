# 必須列を確認する: 指定された変数名が `data` に存在するか検査する。
# Check required columns: verifies that requested variables exist in `data`.
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

# entry age が baseline confounders に含まれていないか確認する: attained age を時間軸とする際の重複モデリングと反実仮想予測時の外挿（age_temp < age_at_entry）を警告する。
# Warn if entry age is in baseline confounders: prevents redundant modelling on the attained-age time scale and extrapolation when predicting counterfactual hazards at age_temp < age_at_entry.
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

# 二値曝露のレベルを決める: 参照群と曝露群の値を確定する。
# Resolve binary exposure levels: determines reference and exposed values.
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
