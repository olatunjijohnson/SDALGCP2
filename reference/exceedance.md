# Exceedance probabilities P(risk \> threshold)

Exceedance probabilities P(risk \> threshold)

## Usage

``` r
exceedance(object, thresholds, which = c("adjusted_rr", "relative_risk"))
```

## Arguments

- object:

  an `"SDALGCP2_pred"` object from
  [`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md).

- thresholds:

  numeric vector of thresholds.

- which:

  which quantity: `"adjusted_rr"` (the covariate-adjusted relative risk
  \\\exp(S)\\, default) or `"relative_risk"` (the relative risk
  \\\exp(d'\beta + S)\\).

## Value

a matrix of exceedance probabilities (locations x thresholds).
