# Map a spatio-temporal prediction for one time

Maps a chosen quantity (`"relative_risk"`, `"adjusted_rr"`,
`"relative_risk_se"`, `"adjusted_rr_se"` or `"exceedance"`) for a
selected time slice of a spatio-temporal prediction.

## Usage

``` r
# S3 method for class 'SDALGCP2_ST_pred'
plot(
  x,
  time = x$times[1],
  what = c("relative_risk", "adjusted_rr", "relative_risk_se", "adjusted_rr_se",
    "exceedance"),
  threshold = 1,
  which = c("adjusted_rr", "relative_risk"),
  ...
)
```

## Arguments

- x:

  an `"SDALGCP2_ST_pred"` object from
  [`predict()`](https://rdrr.io/r/stats/predict.html) on an
  `"SDALGCP2_ST"` fit.

- time:

  the time to map (one of the fitted `times`); defaults to the first.
  Use `NULL` to facet all times.

- what:

  one of `"relative_risk"`, `"adjusted_rr"`, `"relative_risk_se"`,
  `"adjusted_rr_se"`, `"exceedance"`.

- threshold:

  threshold for `what = "exceedance"`.

- which:

  for exceedance: `"adjusted_rr"` (default) or `"relative_risk"`.

- ...:

  unused.

## Value

a `ggplot` object.
