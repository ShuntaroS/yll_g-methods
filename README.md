# YLL g-formula prototype

このフォルダは、YLL (years of life lost) を g-formula で推定する R 実装の試作です。

## 現在のファイル

- `data.rds`: toy data
- `yll_g_formula.R`: 推定関数
- `validation/validate_yll_algorithm.R`: 既知モデルと null effect を使った検証スクリプト
- `validation/validate_refactor_equivalence.R`: リファクタ前後の一致検証スクリプト

## このリポジトリの toy data

`data.rds` には次の変数が入っています。

- `id`: 個体 ID
- `period`: 観察開始からイベントまたは打ち切りまでの追跡時間
- `event`: イベント指標
- `smoke_binary`: 喫煙曝露
- `hypertension`: 高血圧曝露
- `sex`: 性別
- `education_years`: 学歴年数
- `bmi`: BMI
- `age`: 観察開始年齢

この toy data は、

- `smoke_binary` と `hypertension` の両方がイベントに影響する
- `sex`, `education_years`, `bmi`, `age` が交絡因子として振る舞う
- `age` は整数

という設定で作っています。

## 現時点で確認したこと

- 生存確率の時点合わせが 1 区間ずれていたため、LE と YLL が系統的にずれる問題を修正しました。
- `validation/validate_yll_algorithm.R` で、既知の離散時間ハザードモデルに対して推定値が真値に近いことを確認できます。
- 同じスクリプトで、曝露効果が 0 のときに YLL がほぼ 0 になることも確認できます。

## 既知の注意点

- attained age を時間軸にして YLL を出すとき、`age_at_entry_var` が個体ごとに異なる場合は left truncation が問題になります。
- 現在の実装はこの点をまだ厳密には扱っていないため、`age_at_entry_var` がばらつくデータに対する推定は暫定的なものとして扱ってください。

## 現時点の課題

- `age_at_entry_var` が個体ごとに異なるときの `time zero` と交絡因子測定時点の不一致をまだ扱っていません。
- 上記に関連して、可変 entry age のデータでは left truncation / delayed entry による選択バイアスが残りえます。
- conditional survival probability を使った YLL 推定の解釈は、固定 entry age の設定を主な対象としており、可変 entry age では limitation が残ります。
- binary exposure については、`ATT` / `ATC` / `ATE` に加えて deterministic / stochastic な population intervention を指定できます。
- 現在の実装は binary exposure を前提としており、多値曝露や連続曝露の介入指定 API は未整備です。

## 今後の拡張候補

- 多値曝露や連続曝露にも拡張できる intervention API にする。
- 曝露変数、アウトカム、交絡因子の指定をより一般化し、再利用しやすい関数インターフェースにする。
- 最終的に R package 化し、GitHub 上で公開・保守できる形にする。

## 使い方

### 事前準備

```r
source("yll_g_formula.R")
d <- readRDS("data.rds")
```

検証スクリプト:

```sh
Rscript validation/validate_yll_algorithm.R
Rscript validation/validate_refactor_equivalence.R
```

### 推定対象の考え方

- `ATE`: 集団全体で平均した効果
- `ATT`: もともと曝露だった人に限定して平均した効果
- `ATC`: もともと非曝露だった人に限定して平均した効果
- `interventional effect`: 2つの介入シナリオを比較した効果

### `target_population` とは

`target_population` は、「最終的にどの集団で平均をとるか」を指定する引数です。

- `target_population = NULL`
  データ全体で平均します。
- `target_population = function(data) data$sex == "Male"`
  男性だけで平均します。
- `target_population = function(data) data$hypertension == "Yes"`
  もともと高血圧だった人だけで平均します。

計算の流れは、

1. まず与えられたデータ全体でハザードモデルを当てる
2. 全員について反実仮想 survival / LE / YLL を作る
3. 最後に `target_population` で指定した集団に絞って平均する

という順です。つまり、`target_population` は「モデルを当てる前にデータを捨てる」引数ではなく、「反実仮想値を作った後にどの集団で平均するか」を決める引数です。

### toy data で使う変数

以下の例では、toy data の次の変数を使っています。

- 曝露: `hypertension` または `smoke_binary`
- アウトカム: `event`
- 交絡因子: `sex`, `education_years`, `bmi`, `age` と、必要に応じてもう一方の曝露変数

### 結果オブジェクトの見方

推定結果は list として返ります。主に使うのは次の 3 つです。

- `res$summary`
  主要な結果表です。基本的には `starting_age`, `yll`, `ci_low`, `ci_high`, `ci_method` を見れば十分です。
- `res$detailed_results`
  信頼区間の詳細や標準誤差も含む表です。
- `res$meta`
  `B`, `seed`, `conf_level` などの設定情報です。

`method` は `"normal"` または `"percentile"` のどちらか 1 つだけを指定します。デフォルトは `"normal"` です。`res$summary$ci_method` を見れば、どちらの方法で計算した区間か分かります。

### binary stochastic intervention の関数の関係

`estimate_yll_gformula_binary_stochastic_vs_natural()` は、自然経過と比較したいときのユーザー向けラッパーです。実際の処理は次の順で進みます。

```text
estimate_yll_gformula_binary_stochastic_vs_natural()
  -> estimate_yll_gformula_binary_stochastic()
    -> yll_make_binary_stochastic_intervention()
    -> estimate_yll_gformula_intervention()
      -> estimate_yll_gformula_engine_single()
```

- `estimate_yll_gformula_binary_stochastic_vs_natural()`
  自然経過を参照シナリオとして固定します。
- `estimate_yll_gformula_binary_stochastic()`
  2つの stochastic intervention を作って汎用介入関数に渡します。
- `yll_make_binary_stochastic_intervention()`
  各個体について、もとの曝露状態に応じた `P(A*=1)` を返す関数を作ります。
- `estimate_yll_gformula_engine_single()`
  ハザードモデルの当てはめ、反実仮想生存曲線の作成、LE と YLL の計算を行うコア関数です。

`estimate_yll_gformula_engine_single()` の中では、まず各個体について「常に非曝露」の生存曲線と「常に曝露」の生存曲線を作り、その後に介入で指定した `P(A*=1)` で混合して、介入下の平均生存曲線を作っています。

### ATE の例

```r
res_ate <- estimate_yll_gformula_ate(
  data = d,
  B = 200,
  method = "normal",
  show_progress = TRUE,
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
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)
```

この例は、集団全体を対象に「全員非高血圧」対「全員高血圧」の YLL を推定しています。

### ATT の例

```r
res_att <- estimate_yll_gformula_att(
  data = d,
  B = 200,
  method = "normal",
  show_progress = TRUE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "smoke_binary",
  reference_level = "Never",
  exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "hypertension", "age"),
  use_future = FALSE
)
```

この例は、もともと喫煙していた人だけを対象にして、「その人たちが非喫煙だった場合」と「喫煙していた場合」の YLL を比較しています。

### ATC の例

```r
res_atc <- estimate_yll_gformula_atc(
  data = d,
  B = 200,
  method = "normal",
  show_progress = TRUE,
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
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)
```

この例は、もともと高血圧でなかった人だけを対象にして、「その人たちがそのまま高血圧でなかった場合」と「高血圧になった場合」の YLL を比較しています。

### 汎用 interventional effect の例

```r
res_int <- estimate_yll_gformula_intervention(
  data = d,
  intervention_reference = "No",
  intervention_exposed = "Yes",
  target_population = NULL,
  B = 200,
  method = "normal",
  show_progress = TRUE,
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
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)
```

この `estimate_yll_gformula_intervention()` は、2つの介入シナリオの YLL を比較する汎用関数です。上の例では `intervention_reference = "No"`、`intervention_exposed = "Yes"` なので、「全員非高血圧」対「全員高血圧」の YLL を推定しており、実質的には ATE と同じ対象を推定しています。

### 曝露群の一部が非曝露になる interventional effect

```r
res_decrease <- estimate_yll_gformula_binary_stochastic_vs_natural(
  data = d,
  prob_exposed_if_unexposed = 0,
  prob_exposed_if_exposed = 0.7,
  target_population = NULL,
  B = 200,
  show_progress = TRUE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "smoke_binary",
  reference_level = "Never",
  exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "hypertension", "age"),
  use_future = FALSE
)
```

この例では、「もともと喫煙していた人の 30% が禁煙する、もともと喫煙していなかった人はそのまま喫煙しない」という介入を、自然経過と比較しています。

### 非曝露群の一部が曝露になる interventional effect

```r
res_increase <- estimate_yll_gformula_binary_stochastic_vs_natural(
  data = d,
  prob_exposed_if_unexposed = 0.2,
  prob_exposed_if_exposed = 1,
  target_population = NULL,
  B = 200,
  show_progress = TRUE,
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
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary", "age"),
  use_future = FALSE
)
```

この例では、「もともと非高血圧だった人の 20% が高血圧になる、もともと高血圧だった人はそのまま高血圧」という介入を、自然経過と比較しています。

### 自分のデータに使うとき

自分のデータでは、少なくとも次を明示的に指定してください。

- `time_var`
- `event_var`
- `exposure_var`
- `age_at_entry_var`
- `confounders_baseline`

公開関数は特定の列名を前提にしない形にしてあります。
