library(tidyverse)
library(Epi)      # Ns()
library(survival) # Surv(), survSplit()
library(future)
library(future.apply)
library(progressr)



yll_make_cut_points <- function(data, time_var, by = 1L) {
  # survSplit の cut は 0 を含めないベクトルを渡す想定なので，ここでは全体を返す
  max_t <- ceiling(max(data[[time_var]], na.rm = TRUE))
  seq.int(from = 0L, to = max_t, by = by)
}

yll_survsplit <- function(data, time_var, event_var, cut_points, episode = "tgroup") {
  # Surv(time, event) ~ . を survSplit に渡す
  surv_formula <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ."))
  survSplit(
    formula  = surv_formula,
    data     = data,
    cut      = cut_points[-1],
    episode  = episode
  )
}

yll_knot_quantiles <- function(x, probs = c(0, 0.25, 0.50, 0.75, 1)) {
  as.numeric(quantile(x, probs = probs, na.rm = TRUE))
}

yll_rhs_confounders <- function(confounders_baseline) {
  if (is.null(confounders_baseline) || length(confounders_baseline) == 0) {
    "1"
  } else {
    paste(confounders_baseline, collapse = " + ")
  }
}

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

yll_warn_if_varying_entry_age <- function(data, age_at_entry_var) {
  if (dplyr::n_distinct(data[[age_at_entry_var]]) > 1L) {
    warning(
      "`",
      age_at_entry_var,
      "` varies across individuals. Estimating YLL on the attained-age scale ",
      "with varying entry ages can be biased by left truncation; treat the ",
      "current result as provisional until that is handled explicitly.",
      call. = FALSE
    )
  }
}

yll_fit_discrete_hazard_glm <- function(
    df_long,
    event_var,
    exposure_var_internal = "expo",
    age_temp_var = "age_temp",
    confounders_baseline = NULL,
    knot_p
) {
  rhs <- yll_rhs_confounders(confounders_baseline)
  
  # 交互作用込み
  form <- as.formula(
    paste0(
      event_var,
      " ~ ",
      exposure_var_internal,
      " + Epi::Ns(", age_temp_var, ", knots = knot_p) + ",
      exposure_var_internal, " * Epi::Ns(", age_temp_var, ", knots = knot_p) + ",
      rhs
    )
  )
  
  glm(formula = form, family = binomial(), data = df_long)
}

yll_expand_counterfactual_data <- function(
    data,
    id_var,
    exposure_var,
    age_end_inclusive,
    exposure_var_internal = "expo"
) {
  # 各個体について age_temp = 0..age_end_inclusive を作成
  # uncount を使うため tidyr 前提
  out <- data |>
    uncount(weights = age_end_inclusive + 1L, .remove = FALSE) |>
    group_by(.data[[id_var]]) |>
    mutate(age_temp = row_number() - 1L) |>
    ungroup() |>
    rename(expo_original = all_of(exposure_var))
  
  # モデル用 exposure 列は別名で持つ
  out[[exposure_var_internal]] <- out$expo_original
  out
}

yll_predict_hazard <- function(fit, newdata) {
  as.numeric(predict(fit, newdata = newdata, type = "response"))
}

yll_select_estimand_population <- function(df0, df1, estimand, reference_level, exposed_level) {
  if (identical(estimand, "ATT")) {
    df0 <- filter(df0, expo_original == exposed_level)
    df1 <- filter(df1, expo_original == exposed_level)
  } else if (identical(estimand, "ATC")) {
    df0 <- filter(df0, expo_original == reference_level)
    df1 <- filter(df1, expo_original == reference_level)
  } else if (!identical(estimand, "ATE")) {
    stop("estimand must be one of 'ATT', 'ATC', or 'ATE'.")
  }
  list(df0 = df0, df1 = df1)
}

yll_mean_survival_by_age <- function(df, id_var, hazard_var, surv_name = "surv") {
  # Survival is indexed at the start of each age interval.
  df |>
    arrange(.data[[id_var]], age_temp) |>
    group_by(.data[[id_var]]) |>
    mutate(
      .surv_after_interval = cumprod(1 - .data[[hazard_var]]),
      .surv_at_age = dplyr::lag(.surv_after_interval, default = 1)
    ) |>
    ungroup() |>
    group_by(age_temp) |>
    summarise(!!surv_name := mean(.surv_at_age), .groups = "drop")
}

yll_conditional_survival_from_age <- function(g_surv_data, a_start, surv0_name = "surv0", surv1_name = "surv1") {
  s0 <- g_surv_data[[surv0_name]][g_surv_data$age_temp == a_start]
  s1 <- g_surv_data[[surv1_name]][g_surv_data$age_temp == a_start]
  
  if (length(s0) == 0 || length(s1) == 0 || is.na(s0) || is.na(s1) || s0 <= 0 || s1 <= 0) {
    return(NULL)
  }
  
  g_surv_data |>
    filter(age_temp >= a_start) |>
    mutate(
      surv_cond_g0 = .data[[surv0_name]] / s0,
      surv_cond_g1 = .data[[surv1_name]] / s1
    )
}

yll_integrate_le <- function(age, surv) {
  if (length(age) <= 1) {
    return(0)
  }
  
  ord <- order(age)
  age <- age[ord]
  surv <- surv[ord]
  
  sum(diff(age) * head(surv, -1))
}

etimate_yll_gformula_single <- function(
    data,
    id_var,
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start,
    age_end,
    age_interval,
    confounders_baseline = NULL,
    estimand = c("ATT", "ATC", "ATE"),
    cut_by = 1L
) {
  estimand <- match.arg(estimand)
  
  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  
  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level
  
  # survSplit 用のカット
  cut_points <- yll_make_cut_points(data, time_var = time_var, by = cut_by)
  
  # 離散時間ロング化
  df_long <- yll_survsplit(
    data       = data,
    time_var   = time_var,
    event_var  = event_var,
    cut_points = cut_points
  ) |>
    mutate(
      age_temp = .data[[age_at_entry_var]] + tstart,
      expo     = .data[[exposure_var]]
    )
  
  knot_p <- yll_knot_quantiles(df_long$age_temp)
  
  fit <- yll_fit_discrete_hazard_glm(
    df_long              = df_long,
    event_var            = event_var,
    exposure_var_internal = "expo",
    age_temp_var         = "age_temp",
    confounders_baseline = confounders_baseline,
    knot_p               = knot_p
  )
  
  # 反実仮想データ（年齢 0..age_end を作る）
  age_end_inclusive <- as.integer(age_end)
  cf_base <- yll_expand_counterfactual_data(
    data              = data,
    id_var            = id_var,
    exposure_var      = exposure_var,
    age_end_inclusive = age_end_inclusive,
    exposure_var_internal = "expo"
  )
  
  cf0 <- cf_base
  cf1 <- cf_base
  cf0$expo <- reference_level
  cf1$expo <- exposed_level
  
  cf0$hazard0 <- yll_predict_hazard(fit, cf0)
  cf1$hazard1 <- yll_predict_hazard(fit, cf1)
  
  # 対象集団（ATT / ATC / ATE）
  pop <- yll_select_estimand_population(cf0, cf1, estimand, reference_level, exposed_level)
  cf0 <- pop$df0
  cf1 <- pop$df1
  
  # 生存曲線（年齢ごと平均）
  s0 <- yll_mean_survival_by_age(cf0, id_var = id_var, hazard_var = "hazard0", surv_name = "surv0")
  s1 <- yll_mean_survival_by_age(cf1, id_var = id_var, hazard_var = "hazard1", surv_name = "surv1")
  
  g_surv <- inner_join(s0, s1, by = "age_temp")
  
  age_list <- seq(from = age_start, to = age_end, by = age_interval)
  
  map_dfr(age_list, function(a_start) {
    df_cond <- yll_conditional_survival_from_age(g_surv, a_start, "surv0", "surv1")
    if (is.null(df_cond)) {
      return(tibble(age_start = a_start, yll = NA_real_, le_m0 = NA_real_, le_m1 = NA_real_))
    }
    
    le_m0 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g0)
    le_m1 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g1)
    
    tibble(
      age_start = a_start,
      yll  = le_m0 - le_m1,
      le_m0 = le_m0,
      le_m1 = le_m1
    )
  })
}

# ------------------------------------------------------------------------------
# Bootstrap
# ------------------------------------------------------------------------------

yll_bootstrap_ids <- function(ids) {
  n <- length(ids)
  sample(ids, size = n, replace = TRUE)
}

yll_make_boot_data <- function(data, id_var, sampled_ids) {
  # 同一 id が複数回サンプルされた場合に備えて一意な boot_id を振る
  id_map <- tibble(
    !!id_var := sampled_ids,
    boot_id  = paste0(sampled_ids, "_rep", ave(sampled_ids, sampled_ids, FUN = seq_along))
  )
  
  left_join(id_map, data, by = id_var) |>
    mutate(!!id_var := boot_id) |>
    select(-boot_id)
}

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

estimate_yll_gformula <- function(
    data,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    method = c("percentile", "normal"),
    show_progress = TRUE,
    id_var = "id",
    time_var = "period",
    event_var = "dflag",
    exposure_var = "smoke_binary",
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var = "age_in",
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = c("male", "grad_high", "age_in"),
    estimand = c("ATT", "ATC", "ATE"),
    # 並列は呼び出し側で future::plan() を設定する前提
    use_future = TRUE
) {
  estimand <- match.arg(estimand)
  method <- unique(method)
  
  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  yll_warn_if_varying_entry_age(data, age_at_entry_var)
  
  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level
  
  set.seed(seed)
  
  ids <- unique(data[[id_var]])
  
  point_est <- etimate_yll_gformula_single(
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
    estimand = estimand
  )
  
  if (B < 1L) {
    if (length(method) > 0) {
      warning("B = 0 so confidence intervals were not computed.", call. = FALSE)
    }
    
    return(list(
      point   = point_est,
      boot    = tibble(),
      summary = arrange(point_est, age_start),
      meta    = list(B = B, seed = seed, conf_level = conf_level, method = method, estimand = estimand)
    ))
  }
  
  one_boot <- function(b) {
    sampled_ids <- yll_bootstrap_ids(ids)
    boot_data <- yll_make_boot_data(data, id_var, sampled_ids)
    
    etimate_yll_gformula_single(
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
      estimand = estimand
    ) |>
      mutate(b = b)
  }
  
  if (use_future) {
    if (show_progress) {
      handlers(global = TRUE)
      handlers("txtprogressbar")
      boot_list <- with_progress({
        p <- progressor(steps = B)
        future_lapply(seq_len(B), function(b) { p(); one_boot(b) }, future.seed = TRUE)
      })
    } else {
      boot_list <- future_lapply(seq_len(B), one_boot, future.seed = TRUE)
    }
  } else {
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)
      boot_list <- lapply(seq_len(B), function(b) { setTxtProgressBar(pb, b); one_boot(b) })
    } else {
      boot_list <- lapply(seq_len(B), one_boot)
    }
  }
  
  boot_df <- bind_rows(boot_list)
  
  summary <- point_est
  
  if ("percentile" %in% method) {
    perc <- yll_ci_percentile(boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_perc"), -age_start)
    
    summary <- left_join(summary, perc, by = "age_start")
  }
  
  if ("normal" %in% method) {
    norm <- yll_ci_normal(point_est, boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_norm"), -age_start)
    
    summary <- left_join(summary, norm, by = "age_start")
  }
  
  summary <- arrange(summary, age_start)
  
  list(
    point   = point_est,
    boot    = boot_df,
    summary = summary,
    meta    = list(B = B, seed = seed, conf_level = conf_level, method = method, estimand = estimand)
  )
}
