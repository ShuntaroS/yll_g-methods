source("yll_g_formula.R")
d <- readRDS("data.rds")

res_ate <- estimate_yll_gformula_ate(
  data = d,
  B = 50,
  method = "normal",
  show_progress = TRUE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)


res_ate


res_att <- estimate_yll_gformula_att(
  data = d,
  B = 50,
  method = "normal",
  show_progress = TRUE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)

res_att


res_atc <- estimate_yll_gformula_atc(
  data = d,
  B = 50,
  method = "normal",
  show_progress = TRUE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)



