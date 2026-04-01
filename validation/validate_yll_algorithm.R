source("yll_g_formula.R")

run_known_model_check <- function(seed = 42, n = 40000) {
  set.seed(seed)
  
  max_follow <- 5
  age0 <- rep(50, n)
  expo <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("Never", "Current/Ever"))
  
  alpha <- -4
  beta_a <- 0.6
  beta_age <- 0.25
  
  period <- rep(max_follow, n)
  dflag <- integer(n)
  
  for (i in seq_len(n)) {
    for (j in 0:(max_follow - 1)) {
      hazard_ij <- plogis(alpha + beta_a * (expo[i] == "Current/Ever") + beta_age * j)
      if (runif(1) < hazard_ij) {
        period[i] <- j + 1
        dflag[i] <- 1L
        break
      }
    }
  }
  
  sim_data <- data.frame(
    id = seq_len(n),
    period = period,
    dflag = dflag,
    smoke_binary = expo,
    age = age0
  )
  
  est <- estimate_yll_gformula(
    data = sim_data,
    B = 0,
    seed = seed,
    method = character(0),
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "dflag",
    exposure_var = "smoke_binary",
    reference_level = "Never",
    exposed_level = "Current/Ever",
    age_at_entry_var = "age",
    age_start = 50,
    age_end = 55,
    age_interval = 5,
    confounders_baseline = NULL,
    estimand = "ATE",
    use_future = FALSE
  )$point |>
    filter(age_start == 50)
  
  haz0 <- sapply(0:(max_follow - 1), function(j) plogis(alpha + beta_age * j))
  haz1 <- sapply(0:(max_follow - 1), function(j) plogis(alpha + beta_a + beta_age * j))
  surv0 <- cumprod(c(1, 1 - haz0))
  surv1 <- cumprod(c(1, 1 - haz1))
  
  truth <- tibble(
    truth_le_m0 = sum(surv0[1:max_follow]),
    truth_le_m1 = sum(surv1[1:max_follow]),
    truth_yll = truth_le_m0 - truth_le_m1
  )
  
  bind_cols(est, truth) |>
    mutate(
      abs_err_yll = abs(yll - truth_yll),
      abs_err_le_m0 = abs(le_m0 - truth_le_m0),
      abs_err_le_m1 = abs(le_m1 - truth_le_m1)
    )
}

run_null_effect_check <- function(seed = 7, n = 40000) {
  set.seed(seed)
  
  max_follow <- 5
  age0 <- rep(50, n)
  expo <- factor(rbinom(n, 1, 0.5), levels = c(0, 1), labels = c("Never", "Current/Ever"))
  
  alpha <- -4
  beta_age <- 0.25
  
  period <- rep(max_follow, n)
  dflag <- integer(n)
  
  for (i in seq_len(n)) {
    for (j in 0:(max_follow - 1)) {
      hazard_ij <- plogis(alpha + beta_age * j)
      if (runif(1) < hazard_ij) {
        period[i] <- j + 1
        dflag[i] <- 1L
        break
      }
    }
  }
  
  sim_data <- data.frame(
    id = seq_len(n),
    period = period,
    dflag = dflag,
    smoke_binary = expo,
    age = age0
  )
  
  estimate_yll_gformula(
    data = sim_data,
    B = 0,
    seed = seed,
    method = character(0),
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "dflag",
    exposure_var = "smoke_binary",
    reference_level = "Never",
    exposed_level = "Current/Ever",
    age_at_entry_var = "age",
    age_start = 50,
    age_end = 55,
    age_interval = 5,
    confounders_baseline = NULL,
    estimand = "ATE",
    use_future = FALSE
  )$point |>
    filter(age_start == 50) |>
    transmute(abs_yll = abs(yll))
}

known_model <- run_known_model_check()
null_effect <- run_null_effect_check()

cat("Known-model validation:\n")
print(known_model)
cat("\nNull-effect validation:\n")
print(null_effect)

stopifnot(known_model$abs_err_yll < 0.05)
stopifnot(known_model$abs_err_le_m0 < 0.05)
stopifnot(known_model$abs_err_le_m1 < 0.05)
stopifnot(null_effect$abs_yll < 0.05)

cat("\nAll validation checks passed.\n")
