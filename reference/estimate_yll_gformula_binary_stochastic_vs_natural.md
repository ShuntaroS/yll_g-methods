# Compare a binary stochastic intervention against the natural course

Convenience wrapper around
[`estimate_yll_gformula_binary_stochastic()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_binary_stochastic.md)
in which the *reference* scenario is fixed to the natural course (every
subject keeps their observed exposure).

## Usage

``` r
estimate_yll_gformula_binary_stochastic_vs_natural(
  data,
  prob_exposed_if_unexposed,
  prob_exposed_if_exposed,
  target_population = NULL,
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
  integration = c("left_rectangle", "trapezoidal"),
  use_future = TRUE
)
```

## Arguments

- data:

  A `data.frame` with one row per individual, containing the variables
  specified by `id_var`, `time_var`, `event_var`, `exposure_var`,
  `age_at_entry_var`, and any `confounders_baseline`.

- prob_exposed_if_unexposed, prob_exposed_if_exposed:

  Probabilities defining the comparison intervention.

- target_population:

  Either `NULL` (use the full population) or a function `f(data)` that
  returns a logical vector marking the subjects to average over.
  Subjects can be selected on `data$expo_original`.

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

Same list structure as
[`estimate_yll_gformula()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula.md).
