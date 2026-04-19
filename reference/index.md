# Package index

## Estimation

Top-level functions for estimating YLL.

- [`estimate_yll_gformula()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula.md)
  : Estimate YLL via the g-formula (ATE / ATT / ATC)
- [`estimate_yll_gformula_ate()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_estimand_wrappers.md)
  [`estimate_yll_gformula_att()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_estimand_wrappers.md)
  [`estimate_yll_gformula_atc()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_estimand_wrappers.md)
  : Estimand-specific wrappers
- [`estimate_yll_gformula_intervention()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_intervention.md)
  : Estimate YLL under arbitrary interventions and target population
- [`estimate_yll_gformula_binary_stochastic()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_binary_stochastic.md)
  : Compare two binary stochastic interventions
- [`estimate_yll_gformula_binary_stochastic_vs_natural()`](https://shuntaros.github.io/yll_g-methods/reference/estimate_yll_gformula_binary_stochastic_vs_natural.md)
  : Compare a binary stochastic intervention against the natural course

## Intervention helpers

- [`yll_make_binary_stochastic_intervention()`](https://shuntaros.github.io/yll_g-methods/reference/yll_make_binary_stochastic_intervention.md)
  : Build a binary stochastic intervention

## Plotting

- [`plot_yll()`](https://shuntaros.github.io/yll_g-methods/reference/plot_yll.md)
  : Plot YLL across starting ages
- [`plot_marginal_survival()`](https://shuntaros.github.io/yll_g-methods/reference/plot_marginal_survival.md)
  : Plot marginal survival curves for the two intervention arms
- [`plot_conditional_survival()`](https://shuntaros.github.io/yll_g-methods/reference/plot_conditional_survival.md)
  : Plot the conditional survival probability from a starting age

## Data

- [`yll_toy`](https://shuntaros.github.io/yll_g-methods/reference/yll_toy.md)
  : Synthetic survival dataset for YLL examples
