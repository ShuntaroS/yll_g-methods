# Draw one bootstrap sample of *individual* IDs (not rows). Resampling at the
# subject level is essential here because the modelling pipeline treats each
# subject as a unit and would otherwise underestimate uncertainty by treating
# correlated person-period rows as independent.
#' @noRd
yll_bootstrap_ids <- function(ids) {
  n <- length(ids)
  sample(ids, size = n, replace = TRUE)
}

# Materialise a bootstrap dataset from the resampled ID vector.
#
# Naively `left_join`ing on the original ID column would give duplicated rows
# that the engine would (incorrectly) treat as the same subject. To keep each
# resampled copy independent we mint a fresh "_repK" suffix per duplicate,
# then assign that synthetic ID back to the original column so downstream
# code is none the wiser.
#' @noRd
yll_make_boot_data <- function(data, id_var, sampled_ids) {
  id_map <- tibble(
    !!id_var := sampled_ids,
    boot_id  = paste0(sampled_ids, "_rep", ave(sampled_ids, sampled_ids, FUN = seq_along))
  )

  left_join(id_map, data, by = id_var) |>
    mutate(!!id_var := .data[["boot_id"]]) |>
    select(-"boot_id")
}

# Percentile bootstrap CIs for YLL and both arm-specific LEs.
#
# `boot_df` is a long tibble with one row per (bootstrap iteration,
# starting age). For each starting age we take the requested lower/upper
# quantiles of the bootstrap distribution.
#' @noRd
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

# Wald-type ("normal-approximation") bootstrap CIs.
#
# Standard errors come from the bootstrap distribution; the interval is
# point_est ± z_{1-alpha/2} * SE_boot. This is symmetric around the point
# estimate and is the right choice when the sampling distribution is roughly
# normal — fast to compute and well-behaved at moderate B.
#' @noRd
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

# Rename the internal column names (`yll_lwr_perc`, `se_yll_norm`, etc.) into
# more readable user-facing labels for the result object's `detailed_results`
# table, and assemble the trimmed `summary` table whose `ci_low` / `ci_high`
# columns reflect the user's chosen `method`.
#
# Two tables are returned:
#   * `detailed_results` — every CI / SE column under both methods; useful for
#     downstream analysis and sanity checks.
#   * `summary`          — one row per starting age with point estimate, the
#     requested CI, and the CI method tag. This is what most users actually
#     read.
#' @noRd
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

  # Default to NA bounds; the branches below fill them in with whichever CI
  # method the user requested.
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
        ci_low = .data[["ci_low_percentile"]],
        ci_high = .data[["ci_high_percentile"]],
        ci_method = method
      )
  } else if (identical(method, "normal") &&
             all(c("ci_low_normal", "ci_high_normal") %in% names(detailed_results))) {
    main_results <- detailed_results |>
      transmute(
        starting_age,
        yll,
        ci_low = .data[["ci_low_normal"]],
        ci_high = .data[["ci_high_normal"]],
        ci_method = method
      )
  }

  list(
    detailed_results = detailed_results,
    summary = main_results
  )
}

# Assemble the final result list returned to the user.
#
# Returning a structured list (instead of dropping a single tibble back) lets
# us carry the run metadata, both summaries, and the marginal survival curves
# needed by the plotting helpers. Plot helpers pick up
# `marginal_survival_point` / `marginal_survival_boot`; downstream analyses
# typically read `summary` and `detailed_results`.
#' @noRd
yll_build_result_object <- function(point_est, boot_df, summary_df, meta, method,
                                    marginal_survival_point = NULL,
                                    marginal_survival_boot = NULL) {
  readable <- yll_make_readable_results(point_est, boot_df, summary_df, method = method)

  list(
    detailed_results = readable$detailed_results,
    summary = readable$summary,
    meta = meta,
    marginal_survival_point = marginal_survival_point,
    marginal_survival_boot  = marginal_survival_boot
  )
}

# Stack the per-iteration marginal survival curves into one long tibble
# (with a `b` column identifying the iteration). NULL-safe and length-safe so
# that bootstrap loops with degenerate iterations don't blow up.
#' @noRd
yll_collect_bootstrap_curves <- function(boot_results) {
  curves <- lapply(boot_results, `[[`, "curves")
  curves <- curves[!vapply(curves, is.null, logical(1))]
  if (length(curves) == 0) return(NULL)
  bind_rows(curves)
}

# Package one bootstrap iteration's outputs into the shape the parent loop
# expects: a list with the YLL tibble (tagged with iteration index `b`) and
# the iteration's marginal survival curves (also tagged with `b`).
#
# The marginal curves come back as an attribute on the YLL tibble so the
# inner engine can return them without changing its return type; here is
# where we promote them into a proper top-level element.
#' @noRd
yll_one_boot_result <- function(yll_tibble, b) {
  curves <- attr(yll_tibble, "marginal_curves")
  if (!is.null(curves)) {
    curves[["b"]] <- b
  }
  list(yll = mutate(yll_tibble, b = b), curves = curves)
}
