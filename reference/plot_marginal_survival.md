# Plot marginal survival curves for the two intervention arms

Draws the population-level marginal survival curves \\S(t)\\ for the two
intervention arms, optionally with pointwise bootstrap confidence bands.

## Usage

``` r
plot_marginal_survival(
  res,
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
