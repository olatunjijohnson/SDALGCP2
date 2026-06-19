# Map a spatio-temporal prediction for one time

Maps a chosen quantity (`"relative_risk"`, `"adjusted_rr"`,
`"relative_risk_se"`, `"adjusted_rr_se"` or `"exceedance"`) for a
selected time slice of a spatio-temporal prediction.

## Usage

``` r
# S3 method for class 'SDALGCP2_ST_pred'
plot(
  x,
  time = attr(x, "times")[1],
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

## Examples

``` r
# \donttest{
data(sdalgcp_data)
times <- 1:3
panel <- do.call(rbind, lapply(times, function(t) {
  d <- sdalgcp_data; d$time <- t
  d$cases <- rpois(nrow(d), d$pop * exp(-6 + 0.6 * d$x1 + 0.1 * (t - 2)))
  d
}))
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = panel, time = "time",
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
pr <- predict(fit)
plot(pr, time = 2)                 # one time slice

plot(pr, time = NULL)              # facet all times

plot(pr, what = "exceedance", threshold = 1.2, time = 3)

# }
```
