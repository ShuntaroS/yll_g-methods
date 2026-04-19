# Synthetic survival dataset for YLL examples

A simulated cohort used to demonstrate the package functions. Each row
is one individual followed for up to 15 years, with a binary smoking
exposure, a binary hypertension covariate, and a
continuous-time-to-event outcome. The data-generating process is encoded
in `data-raw/yll_toy.R`.

## Usage

``` r
yll_toy
```

## Format

A data frame with 5,000 rows and 9 variables:

- id:

  Integer subject identifier (1..n).

- period:

  Numeric follow-up time, in years, from study entry to event or
  administrative censoring at year 15.

- event:

  Integer event indicator: `1` if death, `0` if censored.

- smoke_binary:

  Factor exposure with levels `"Never"` and `"Current/Ever"`.

- hypertension:

  Factor confounder with levels `"No"` and `"Yes"`.

- sex:

  Factor with levels `"Female"` and `"Male"`.

- education_years:

  Integer years of completed education (9-20).

- bmi:

  Numeric body-mass index at entry (16.0-40.0).

- age:

  Integer age in years at study entry (35-69).

## Source

Simulated. See `data-raw/yll_toy.R` for the data-generating script.
