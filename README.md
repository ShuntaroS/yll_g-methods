# YLLgmethods

Years of life lost (YLL) and counterfactual life expectancy under
user-specified interventions on a binary exposure, estimated with the
parametric **g-formula** on the attained-age time scale.

The package supports:

- ATE / ATT / ATC contrasts for a binary exposure.
- Arbitrary deterministic, stochastic, or covariate-dependent
  interventions through a small intervention API.
- Bootstrap confidence intervals (normal-approximation or percentile),
  optionally parallelised with `future`.
- Plotting helpers for marginal and conditional survival curves and for
  YLL across starting ages.

## Installation

The package is currently distributed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ShuntaroS/yll_g-methods")
```

## Quick start

```r
library(YLLgmethods)
data(yll_toy)

res <- estimate_yll_gformula_ate(
  data = yll_toy,
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
  B = 200,
  method = "normal",
  use_future = FALSE
)

res$summary
plot_yll(res)
```

See `vignette("getting-started", package = "YLLgmethods")` for a tour of
the user-facing functions, and
`vignette("interventional-effects", package = "YLLgmethods")` for the
identification framework behind stochastic and dynamic interventions.

## Result object

Every `estimate_yll_gformula_*()` function returns a list with:

| Element                    | Description                                                                  |
| -------------------------- | ---------------------------------------------------------------------------- |
| `summary`                  | One row per starting age: point estimate, CI, CI method.                     |
| `detailed_results`         | Same plus per-arm life expectancy and both percentile and normal CIs.        |
| `meta`                     | Call metadata (`B`, `seed`, `conf_level`, `integration`, `estimand`, …).     |
| `marginal_survival_point`  | Population-level marginal survival curves under each intervention arm.       |
| `marginal_survival_boot`   | Same, one row per bootstrap iteration (used by the plotting helpers).        |

`method` selects between normal-approximation (`"normal"`, default) and
percentile (`"percentile"`) bootstrap CIs; `res$summary$ci_method`
records which one was used.

## API at a glance

Estimands:

- `estimate_yll_gformula_ate()` / `_att()` / `_atc()` — binary-exposure
  contrasts averaged over the chosen target population.
- `estimate_yll_gformula_intervention()` — generic interventional
  contrast accepting either fixed exposure values or a function
  `f(data)` returning subject-specific exposure probabilities.
- `estimate_yll_gformula_binary_stochastic()` and
  `estimate_yll_gformula_binary_stochastic_vs_natural()` — convenience
  wrappers for binary stochastic shifts, the latter using the natural
  exposure distribution as the reference.

Helpers:

- `yll_make_binary_stochastic_intervention()` — build an intervention
  function from `(P(A^*=\text{exposed}\mid A=\text{ref}), P(A^*=\text{exposed}\mid A=\text{exp}))`.
- `plot_yll()`, `plot_marginal_survival()`, `plot_conditional_survival()`
  — `ggplot2`-based plots of the result object.

## Caveats

- The current implementation is for **binary** exposures. Multi-level
  and continuous exposures are not yet supported.
- Variable entry age introduces a left-truncation / delayed-entry
  problem that the current implementation does not address rigorously;
  treat results from datasets with variable `age_at_entry_var` as
  preliminary.

## License

MIT © 2026 Shuntaro Sato.
