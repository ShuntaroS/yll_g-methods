# Generate the `yll_toy` example dataset shipped with the package.
#
# Run with:
#   source("data-raw/yll_toy.R")

set.seed(20260401)

n <- 5000L
max_follow <- 15L

age <- round(rnorm(n, mean = 54.5, sd = 7))
age <- pmin(pmax(age, 35), 69)

sex <- factor(
  rbinom(n, size = 1, prob = 0.52),
  levels = c(0, 1),
  labels = c("Female", "Male")
)

education_years <- round(
  rnorm(
    n,
    mean = 14.2 - 0.3 * (sex == "Male") - 0.05 * (age - 55),
    sd = 2.2
  )
)
education_years <- pmin(pmax(education_years, 9), 20)

bmi_mean <- 23.5 + 0.06 * (age - 55) + 0.9 * (sex == "Male") - 0.12 * (education_years - 14)
bmi <- round(rnorm(n, mean = bmi_mean, sd = 3.2), 1)
bmi <- pmin(pmax(bmi, 16.0), 40.0)

smoke_lp <- -0.35 +
  0.65 * (sex == "Male") -
  0.14 * (education_years - 14) +
  0.03 * (age - 55) +
  0.02 * (bmi - 25)

smoke_binary <- factor(
  rbinom(n, size = 1, prob = plogis(smoke_lp)),
  levels = c(0, 1),
  labels = c("Never", "Current/Ever")
)

hypertension_lp <- -0.8 +
  0.09 * (age - 55) +
  0.55 * (sex == "Male") +
  0.09 * (bmi - 25) -
  0.08 * (education_years - 14) +
  0.30 * (smoke_binary == "Current/Ever")

hypertension <- factor(
  rbinom(n, size = 1, prob = plogis(hypertension_lp)),
  levels = c(0, 1),
  labels = c("No", "Yes")
)

intervals <- 0:(max_follow - 1L)
age_mat <- outer(age, intervals, "+")

hazard_lp <- -5.8 +
  0.095 * (age_mat - 55) +
  0.42 * (smoke_binary == "Current/Ever") +
  0.60 * (hypertension == "Yes") +
  0.16 * (sex == "Male") +
  0.025 * (bmi - 25) -
  0.05 * (education_years - 14) +
  0.12 * (smoke_binary == "Current/Ever") * (hypertension == "Yes")

hazard_mat <- plogis(hazard_lp)
u_mat <- matrix(runif(n * max_follow), nrow = n, ncol = max_follow)
event_mat <- u_mat < hazard_mat

first_event_interval <- apply(event_mat, 1, function(x) {
  event_idx <- which(x)[1]
  if (is.na(event_idx)) {
    max_follow
  } else {
    event_idx
  }
})

event <- as.integer(rowSums(event_mat) > 0)
period <- rep(max_follow, n)

if (any(event == 1L)) {
  period[event == 1L] <- (first_event_interval[event == 1L] - 1) + runif(sum(event == 1L))
}

yll_toy <- data.frame(
  id = seq_len(n),
  period = period,
  event = event,
  smoke_binary = smoke_binary,
  hypertension = hypertension,
  sex = sex,
  education_years = education_years,
  bmi = bmi,
  age = age
)

usethis::use_data(yll_toy, overwrite = TRUE)
