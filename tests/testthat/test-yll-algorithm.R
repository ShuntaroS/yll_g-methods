# Algorithmic validation: the discrete-time hazard g-formula should recover
# the YLL implied by a known data-generating mechanism, and should produce a
# null effect when the exposure has no effect on the hazard.

simulate_known_model <- function(seed, n,
                                 alpha = -4,
                                 beta_a = 0.6,
                                 beta_age = 0.25,
                                 max_follow = 5L) {
  set.seed(seed)
  age0 <- rep(50, n)
  expo <- factor(
    rbinom(n, 1, 0.5),
    levels = c(0, 1),
    labels = c("Never", "Current/Ever")
  )

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

  list(
    data = data.frame(
      id = seq_len(n),
      period = period,
      dflag = dflag,
      smoke_binary = expo,
      age = age0
    ),
    truth = {
      haz0 <- vapply(0:(max_follow - 1), function(j) plogis(alpha + beta_age * j), numeric(1))
      haz1 <- vapply(0:(max_follow - 1), function(j) plogis(alpha + beta_a + beta_age * j), numeric(1))
      surv0 <- cumprod(c(1, 1 - haz0))
      surv1 <- cumprod(c(1, 1 - haz1))
      list(
        le_m0 = sum(surv0[seq_len(max_follow)]),
        le_m1 = sum(surv1[seq_len(max_follow)]),
        yll   = sum(surv0[seq_len(max_follow)]) - sum(surv1[seq_len(max_follow)])
      )
    },
    max_follow = max_follow
  )
}

test_that("g-formula recovers YLL from a known discrete-time model", {
  skip_on_cran()

  sim <- simulate_known_model(seed = 42, n = 40000)

  est <- suppressWarnings(estimate_yll_gformula(
    data = sim$data,
    B = 0,
    seed = 42,
    method = "normal",
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "dflag",
    exposure_var = "smoke_binary",
    reference_level = "Never",
    exposed_level = "Current/Ever",
    age_at_entry_var = "age",
    age_start = 50,
    age_end = 50 + sim$max_follow,
    age_interval = sim$max_follow,
    confounders_baseline = NULL,
    estimand = "ATE",
    use_future = FALSE
  ))$detailed_results

  est <- est[est$starting_age == 50, , drop = FALSE]

  expect_lt(abs(est$yll          - sim$truth$yll),   0.05)
  expect_lt(abs(est$le_reference - sim$truth$le_m0), 0.05)
  expect_lt(abs(est$le_exposed   - sim$truth$le_m1), 0.05)
})

test_that("g-formula returns ~0 YLL when exposure has no effect", {
  skip_on_cran()

  set.seed(7)
  n <- 40000
  max_follow <- 5L
  alpha <- -4
  beta_age <- 0.25

  age0 <- rep(50, n)
  expo <- factor(
    rbinom(n, 1, 0.5),
    levels = c(0, 1),
    labels = c("Never", "Current/Ever")
  )

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

  est <- suppressWarnings(estimate_yll_gformula(
    data = sim_data,
    B = 0,
    seed = 7,
    method = "normal",
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "dflag",
    exposure_var = "smoke_binary",
    reference_level = "Never",
    exposed_level = "Current/Ever",
    age_at_entry_var = "age",
    age_start = 50,
    age_end = 50 + max_follow,
    age_interval = max_follow,
    confounders_baseline = NULL,
    estimand = "ATE",
    use_future = FALSE
  ))$detailed_results

  est <- est[est$starting_age == 50, , drop = FALSE]
  expect_lt(abs(est$yll), 0.05)
})
