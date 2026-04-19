# Quick smoke test: a small run on the bundled toy dataset should produce a
# well-formed result object without errors.

test_that("estimate_yll_gformula returns the expected object on toy data", {
  data(yll_toy, envir = environment())

  res <- suppressWarnings(estimate_yll_gformula(
    data = yll_toy,
    id_var = "id",
    time_var = "period",
    event_var = "event",
    exposure_var = "smoke_binary",
    reference_level = "Never",
    exposed_level = "Current/Ever",
    age_at_entry_var = "age",
    age_start = 60,
    age_end = 65,
    age_interval = 5,
    confounders_baseline = c("sex", "education_years", "bmi"),
    estimand = "ATE",
    B = 0,
    show_progress = FALSE,
    use_future = FALSE
  ))

  expect_named(res, c(
    "detailed_results", "summary", "meta",
    "marginal_survival_point", "marginal_survival_boot"
  ))
  expect_s3_class(res$summary, "data.frame")
  expect_true(all(c("starting_age", "yll", "ci_low", "ci_high", "ci_method") %in%
                    names(res$summary)))
  expect_true(all(is.finite(res$summary$yll)))
})
