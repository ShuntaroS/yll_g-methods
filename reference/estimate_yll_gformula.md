# Estimate YLL via the g-formula (ATE / ATT / ATC)

Fits a discrete-time hazard model on the attained-age time scale, builds
counterfactual survival curves under "all reference" vs. "all exposed"
interventions, and reports Years of Life Lost (YLL) and counterfactual
life expectancy (LE) at user-chosen starting ages, with bootstrap
confidence intervals.

## Usage

``` r
estimate_yll_gformula(
  data,
  B = 1000,
  seed = 1,
  conf_level = 0.95,
  method = c("normal", "percentile"),
  show_progress = TRUE,
  id_var = "id",
  time_var,
  event_var,
  exposure_var,
  reference_level = NULL,
  exposed_level = NULL,
  age_at_entry_var,
  age_start = 50,
  age_end = 100,
  age_interval = 5,
  confounders_baseline = NULL,
  estimand = c("ATT", "ATC", "ATE"),
  integration = c("left_rectangle", "trapezoidal"),
  use_future = TRUE
)
```

## Arguments

- data:

  A `data.frame` with one row per individual, containing the variables
  specified by `id_var`, `time_var`, `event_var`, `exposure_var`,
  `age_at_entry_var`, and any `confounders_baseline`.

- B:

  Number of bootstrap iterations. Set to `0` to skip bootstrapping
  (returns point estimates only).

- seed:

  Integer seed for reproducibility.

- conf_level:

  Confidence level for the intervals (e.g. `0.95`).

- method:

  Either `"normal"` (Wald-type interval using bootstrap SE) or
  `"percentile"` (bootstrap percentile interval).

- show_progress:

  Show a progress bar during bootstrap.

- id_var, time_var, event_var, exposure_var:

  Column names in `data`.

- reference_level, exposed_level:

  Values labelling the two exposure levels. If `NULL` they are inferred
  from the data.

- age_at_entry_var:

  Column name for age at study entry.

- age_start, age_end, age_interval:

  Starting ages for which YLL is reported, given as
  `seq(age_start, age_end, age_interval)`.

- confounders_baseline:

  Character vector of baseline confounder column names. Including the
  entry-age column here is discouraged and triggers a warning, since
  attained age is already modelled via splines.

- estimand:

  One of `"ATT"`, `"ATC"`, or `"ATE"`.

- integration:

  Discretisation rule for the area under the conditional survival curve:
  `"left_rectangle"` (default; matches the discrete-time hazard model)
  or `"trapezoidal"`.

- use_future:

  If `TRUE`, parallelise the bootstrap with
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).
  The caller is expected to set a parallel plan (e.g.
  `future::plan(future::multisession())`).

## Value

A list with elements

- `summary`:

  Point estimates and CI for YLL at each starting age.

- `detailed_results`:

  Same as `summary` plus LE under each arm and all CI / SE columns.

- `meta`:

  Run settings (`B`, `seed`, `conf_level`, `method`, `estimand`,
  `age_start`, `age_end`, `age_interval`, `integration`).

- `marginal_survival_point`:

  Population-marginal counterfactual survival curves (point estimate).

- `marginal_survival_boot`:

  Same curves across bootstrap iterations (used by the plotting
  helpers).

## Details

This function is the main legacy interface. For arbitrary interventions
or target populations, use
[`estimate_yll_gformula_intervention()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_intervention.md).

## See also

[`estimate_yll_gformula_intervention()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_intervention.md),
[`estimate_yll_gformula_binary_stochastic_vs_natural()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_binary_stochastic_vs_natural.md),
[`plot_yll()`](https://shuntaros.github.io/yll_g-methods/reference/plot_yll.md).

## Examples

``` r
# \donttest{
data(yll_toy, package = "YLLgmethods")
res <- estimate_yll_gformula_ate(
  data                 = yll_toy,
  B                    = 0,
  show_progress        = FALSE,
  time_var             = "period",
  event_var            = "event",
  exposure_var         = "hypertension",
  reference_level      = "No",
  exposed_level        = "Yes",
  age_at_entry_var     = "age",
  age_start            = 50,
  age_end              = 80,
  age_interval         = 10,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  use_future           = FALSE
)
#> Warning: B = 0 so confidence intervals were not computed.
res$summary
#> # A tibble: 4 × 5
#>   starting_age   yll ci_low ci_high ci_method
#>          <dbl> <dbl>  <dbl>   <dbl> <chr>    
#> 1           50 2.28      NA      NA normal   
#> 2           60 1.78      NA      NA normal   
#> 3           70 0.822     NA      NA normal   
#> 4           80 0         NA      NA normal   
# }
```
