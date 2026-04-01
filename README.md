# YLL g-formula prototype

このフォルダは、YLL (years of life lost) を g-formula で推定する R 実装の試作です。

## 現在のファイル

- `data.rds`: toy data
- `yll_g_formula.R`: 推定関数
- `validation/validate_yll_algorithm.R`: 既知モデルと null effect を使った検証スクリプト

## 現時点で確認したこと

- 生存確率の時点合わせが 1 区間ずれていたため、LE と YLL が系統的にずれる問題を修正しました。
- `validation/validate_yll_algorithm.R` で、既知の離散時間ハザードモデルに対して推定値が真値に近いことを確認できます。
- 同じスクリプトで、曝露効果が 0 のときに YLL がほぼ 0 になることも確認できます。

## 既知の注意点

- attained age を時間軸にして YLL を出すとき、`age_at_entry_var` が個体ごとに異なる場合は left truncation が問題になります。
- 現在の実装はこの点をまだ厳密には扱っていないため、`age_at_entry_var` がばらつくデータに対する推定は暫定的なものとして扱ってください。

## 使い方

検証スクリプト:

```sh
Rscript validation/validate_yll_algorithm.R
```

toy data での実行例:

```r
source("yll_g_formula.R")
d <- readRDS("data.rds")

res <- estimate_yll_gformula(
  data = d,
  B = 100,
  show_progress = TRUE,
  id_var = "id",
  time_var = "period",
  event_var = "dflag",
  exposure_var = "smoke_binary",
  reference_level = "Never",
  exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  confounders_baseline = "age",
  estimand = "ATE",
  use_future = FALSE
)
```
