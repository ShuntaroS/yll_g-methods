# Getting started with YLLgmethods

`YLLgmethods` estimates **years of life lost (YLL)** under the
parametric g-formula with a discrete-time hazard model. This vignette
walks through the typical workflow on the bundled toy dataset.

``` r
library(YLLgmethods)
data(yll_toy)
str(yll_toy)
#> 'data.frame':    5000 obs. of  9 variables:
#>  $ id             : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ period         : num  15 15 15 15 15 ...
#>  $ event          : int  0 0 1 0 0 0 0 0 0 1 ...
#>  $ smoke_binary   : Factor w/ 2 levels "Never","Current/Ever": 1 1 1 2 2 2 2 1 1 2 ...
#>  $ hypertension   : Factor w/ 2 levels "No","Yes": 1 2 1 1 1 1 1 1 2 2 ...
#>  $ sex            : Factor w/ 2 levels "Female","Male": 2 2 2 1 1 2 1 2 1 1 ...
#>  $ education_years: num  14 15 16 11 12 15 18 14 19 14 ...
#>  $ bmi            : num  22 23 22 24.7 26.6 27.6 26.3 25.6 22.7 18.7 ...
#>  $ age            : num  54 49 64 58 47 61 51 57 43 65 ...
```

## ATE: average treatment effect in the whole population

[`estimate_yll_gformula_ate()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_estimand_wrappers.md)
is a thin wrapper that contrasts “everyone reference” vs. “everyone
exposed” and averages over the entire study population.

``` r
res_ate <- estimate_yll_gformula_ate(
  data = yll_toy,
  B = 50,
  method = "normal",
  show_progress = FALSE,
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
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  use_future = FALSE
)

res_ate$summary
```

The result object always contains:

- `summary` — one row per starting age with point estimate, CI, and CI
  method
- `detailed_results` — point estimates plus per-arm life expectancy and
  both percentile and normal-approximation CIs (when `B > 0`)
- `meta` — call metadata (`estimand`, `B`, `conf_level`, `integration`,
  …)
- `marginal_survival_point` and `marginal_survival_boot` — marginal
  survival curves used by the `plot_*()` helpers

## ATT and ATC

`_att()` averages the contrast over those who were observed exposed,
`_atc()` over those observed at the reference level. Same arguments,
different target population:

``` r
res_att <- estimate_yll_gformula_att(
  data = yll_toy,
  B = 50, method = "normal", show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)

res_atc <- estimate_yll_gformula_atc(
  data = yll_toy,
  B = 50, method = "normal", show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)
```

## Stochastic interventions

For policy-relevant questions where “set everyone to exposed” is
unrealistic, use a **stochastic intervention** that shifts the exposure
distribution. The convenience wrapper
[`estimate_yll_gformula_binary_stochastic_vs_natural()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_binary_stochastic_vs_natural.md)
compares a target intervention to the natural exposure distribution
observed in the data.

The example below contrasts a world in which 20% of the unexposed and
100% of the exposed end up exposed against the natural distribution:

``` r
res_increase <- estimate_yll_gformula_binary_stochastic_vs_natural(
  data = yll_toy,
  prob_exposed_if_unexposed = 0.2,
  prob_exposed_if_exposed   = 1.0,
  target_population = NULL,
  B = 50, show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)
```

See
[`vignette("interventional-effects", package = "YLLgmethods")`](https://shuntaros.github.io/yll_g-methods/articles/interventional-effects.md)
for the underlying theory and more examples.

## Visualisation

Three plotting helpers are included; all return `ggplot2` objects so you
can post-process freely.

``` r
plot_yll(res_ate)
plot_conditional_survival(res_ate, age_start = 60)
plot_marginal_survival(res_ate)
```

## Reproducibility tips

- Set `seed` to make the bootstrap deterministic.
- Set `use_future = TRUE` and configure a
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  to parallelise the bootstrap.
- `integration = "trapezoidal"` produces slightly smoother LE estimates;
  the default `"left_rectangle"` is consistent with the discrete-time
  hazard model.
