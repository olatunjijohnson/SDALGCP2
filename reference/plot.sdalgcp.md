# Map an sdalgcp fit

Predicts and maps a chosen quantity. Works for spatial fits (discrete or
continuous) and spatio-temporal fits (select a `time`).

## Usage

``` r
# S3 method for class 'sdalgcp'
plot(
  x,
  what = c("RR", "ARR", "RR_se", "ARR_se", "exceedance"),
  type = c("discrete", "continuous"),
  time = NULL,
  threshold = 1,
  which = c("ARR", "RR"),
  cellsize = NULL,
  sampler = c("mcmc", "laplace"),
  ...
)
```

## Arguments

- x:

  an `"sdalgcp"` fit.

- what:

  one of `"RR"` (relative risk, default), `"ARR"` (covariate-adjusted
  relative risk), `"RR_se"`, `"ARR_se"` or `"exceedance"`.

- type:

  `"discrete"` (default) or `"continuous"` (spatial fits).

- time:

  for spatio-temporal fits, the time to map (default: first; use `NULL`
  to facet all times).

- threshold:

  threshold for `what = "exceedance"`.

- which:

  for exceedance: `"ARR"` (default) or `"RR"`.

- cellsize:

  grid spacing for `type = "continuous"`.

- sampler:

  `"mcmc"` (default) or `"laplace"`.

- ...:

  passed to the mapping layer.

## Value

a `ggplot` object.
