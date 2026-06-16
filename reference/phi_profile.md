# Profile likelihood and confidence interval for the spatial scale phi

Spline-smoothed profile deviance for `phi`, with the `coverage`-level
confidence interval where the deviance crosses the chi-squared cutoff.
Improves on the loess-based interval of `SDALGCP`.

## Usage

``` r
phi_profile(object, coverage = 0.95, plot = TRUE)
```

## Arguments

- object:

  a fitted `"SDALGCP2"` object.

- coverage:

  confidence level.

- plot:

  logical; draw the deviance curve.

## Value

invisibly, a list with the interval and the smoothed profile; a `ggplot`
is drawn when `plot = TRUE`.
