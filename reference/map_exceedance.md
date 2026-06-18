# Map exceedance probabilities P(risk \> threshold)

Map exceedance probabilities P(risk \> threshold)

## Usage

``` r
map_exceedance(
  x,
  threshold = 1,
  which = c("adjusted_rr", "relative_risk"),
  bound = NULL,
  ...
)
```

## Arguments

- x:

  an `"SDALGCP2_pred"` object.

- threshold:

  a single relative-risk threshold.

- which:

  `"adjusted_rr"` (covariate-adjusted, default) or `"relative_risk"`.

- bound:

  optional `sf` boundary (continuous only).

- ...:

  unused.

## Value

a `ggplot` object.
