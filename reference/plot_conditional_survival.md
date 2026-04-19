# Plot the conditional survival probability from a starting age

Draws the conditional survival curves \\S(y \mid a\_{\text{start}})\\
for the two intervention arms, optionally with a pointwise bootstrap
confidence band.

## Usage

``` r
plot_conditional_survival(
  res,
  age_start,
  conf_band = TRUE,
  reference_label = "Reference",
  exposed_label = "Exposed",
  conf_level = NULL
)
```

## Arguments

- res:

  A result object returned by
  [`estimate_yll_gformula()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula.md)
  (or one of its variants).

- age_start:

  Numeric. The starting age \\a\_{\text{start}}\\ used to condition the
  survival curve.

- conf_band:

  Logical. If `TRUE` (default) and bootstrap curves are available, draws
  a pointwise confidence band around each curve.

- reference_label, exposed_label:

  Character labels for the two intervention arms shown in the legend.

- conf_level:

  Confidence level for the pointwise band. Defaults to the value stored
  in `res$meta$conf_level`, falling back to `0.95`.

## Value

A `ggplot` object.
