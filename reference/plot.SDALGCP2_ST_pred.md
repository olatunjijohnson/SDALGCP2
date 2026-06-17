# Map a spatio-temporal prediction for one time

Maps a chosen quantity (`"RR"`, `"ARR"`, `"RR_se"`, `"ARR_se"` or
`"exceedance"`) for a selected time slice of a spatio-temporal
prediction.

## Usage

``` r
# S3 method for class 'SDALGCP2_ST_pred'
plot(
  x,
  time = x$times[1],
  what = c("RR", "ARR", "RR_se", "ARR_se", "exceedance"),
  threshold = 1,
  which = c("ARR", "RR"),
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

  one of `"RR"`, `"ARR"`, `"RR_se"`, `"ARR_se"`, `"exceedance"`.

- threshold:

  threshold for `what = "exceedance"`.

- which:

  for exceedance: `"ARR"` (default) or `"RR"`.

- ...:

  unused.

## Value

a `ggplot` object.
