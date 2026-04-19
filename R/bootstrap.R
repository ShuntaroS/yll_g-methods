# bootstrap 用 ID を標本抽出する: 個体 ID を復元抽出する。
# Sample bootstrap IDs: resamples individual IDs with replacement.
#' @noRd
yll_bootstrap_ids <- function(ids) {
  n <- length(ids)
  sample(ids, size = n, replace = TRUE)
}

# bootstrap データを作る: 重複抽出された個体に新しい ID を振る。
# Build bootstrap data: assigns fresh IDs to duplicated sampled individuals.
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

# パーセンタイル CI を作る: bootstrap 分布から信頼区間を計算する。
# Compute percentile CIs: derives confidence intervals from bootstrap quantiles.
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

# 正規近似 CI を作る: bootstrap 標準誤差から信頼区間を計算する。
# Compute normal-approximation CIs: derives confidence intervals from bootstrap SEs.
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

# 読みやすい列名に直す: 結果テーブルを用途別に分かりやすい名前へ変換する。
# Make result tables readable: renames output columns to more descriptive names.
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

# 結果オブジェクトを作る: 主結果・詳細表・設定情報を返す。
# Build result object: returns the main table, detailed results, and metadata.
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

# bootstrap 反復の marginal curves を集める: 各反復の周辺生存曲線を縦結合する。
# Collect bootstrap marginal curves: row-binds per-iteration marginal survival tables.
#' @noRd
yll_collect_bootstrap_curves <- function(boot_results) {
  curves <- lapply(boot_results, `[[`, "curves")
  curves <- curves[!vapply(curves, is.null, logical(1))]
  if (length(curves) == 0) return(NULL)
  bind_rows(curves)
}

# 1 回の bootstrap 結果を整える: yll tibble と marginal curves を分離して返す。
# Wrap one bootstrap iteration: returns yll tibble and marginal curves separately.
#' @noRd
yll_one_boot_result <- function(yll_tibble, b) {
  curves <- attr(yll_tibble, "marginal_curves")
  if (!is.null(curves)) {
    curves[["b"]] <- b
  }
  list(yll = mutate(yll_tibble, b = b), curves = curves)
}
