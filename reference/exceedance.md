# Exceedance probabilities P(relative risk \> threshold)

Exceedance probabilities P(relative risk \> threshold)

## Usage

``` r
exceedance(object, thresholds)
```

## Arguments

- object:

  an `"SDALGCP2_pred"` object from
  [`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md).

- thresholds:

  numeric vector of thresholds.

## Value

a matrix of exceedance probabilities (locations x thresholds).
