# Plot YLL across starting ages

Draws estimated years of life lost (YLL) as a function of the starting
age, with an optional confidence band when bootstrap CIs are available
in the result object.

## Usage

``` r
plot_yll(res, conf_band = TRUE)
```

## Arguments

- res:

  A result object returned by
  [`estimate_yll_gformula()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula.md)
  (or one of its variants).

- conf_band:

  Logical. If `TRUE` (default) and `ci_low`/`ci_high` are present in
  `res$summary`, draws a confidence band.

## Value

A `ggplot` object.
