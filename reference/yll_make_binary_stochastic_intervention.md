# Build a binary stochastic intervention

Creates an intervention function suitable for the
`intervention_reference` / `intervention_exposed` arguments of
[`estimate_yll_gformula_intervention()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_intervention.md).
For each subject the returned function reports `P(A* = exposed_level)`
as a function of the subject's originally observed exposure status.

## Usage

``` r
yll_make_binary_stochastic_intervention(
  prob_exposed_if_unexposed,
  prob_exposed_if_exposed,
  reference_level,
  exposed_level
)
```

## Arguments

- prob_exposed_if_unexposed:

  Probability of being assigned the exposed level among subjects whose
  observed exposure was the reference level.

- prob_exposed_if_exposed:

  Probability of being assigned the exposed level among subjects whose
  observed exposure was the exposed level.

- reference_level, exposed_level:

  Values labelling the two exposure levels in the original data.

## Value

A function `f(data)` returning a numeric vector of length `nrow(data)`
whose entries are in `[0, 1]`.
