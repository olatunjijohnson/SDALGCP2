# Map exceedance probabilities P(risk \> threshold)

Map exceedance probabilities P(risk \> threshold)

## Usage

``` r
map_exceedance(x, threshold = 1, which = c("ARR", "RR"), bound = NULL, ...)
```

## Arguments

- x:

  an `"SDALGCP2_pred"` object.

- threshold:

  a single relative-risk threshold.

- which:

  `"ARR"` (covariate-adjusted, default) or `"RR"`.

- bound:

  optional `sf` boundary (continuous only).

- ...:

  unused.

## Value

a `ggplot` object.
