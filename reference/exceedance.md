# Exceedance probabilities P(risk \> threshold)

Exceedance probabilities P(risk \> threshold)

## Usage

``` r
exceedance(object, thresholds, which = c("ARR", "RR"))
```

## Arguments

- object:

  an `"SDALGCP2_pred"` object from
  [`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md).

- thresholds:

  numeric vector of thresholds.

- which:

  which quantity: `"ARR"` (covariate-adjusted relative risk, default) or
  `"RR"` (relative risk).

## Value

a matrix of exceedance probabilities (locations x thresholds).
