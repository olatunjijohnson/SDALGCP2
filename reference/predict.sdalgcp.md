# Predict relative risk from an sdalgcp fit

Returns a prediction object carrying, for every location, the posterior
mean and standard error of the relative risk `RR`
(\\\exp(\eta)=\exp(d'\beta+S)\\) and the covariate-adjusted relative
risk `ARR` (\\\exp(S)\\). Map it with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and get hotspot
probabilities with
[`exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md).

## Usage

``` r
# S3 method for class 'sdalgcp'
predict(
  object,
  type = c("discrete", "continuous"),
  sampler = c("mcmc", "laplace"),
  cellsize = NULL,
  ...
)
```

## Arguments

- object:

  an `"sdalgcp"` fit.

- type:

  `"discrete"` (region level, default) or `"continuous"` (a grid
  surface). Ignored for spatio-temporal fits.

- sampler:

  `"mcmc"` (default) or `"laplace"`.

- cellsize:

  grid spacing for `type = "continuous"`.

- ...:

  passed to the underlying predictor.

## Value

an object of class `"SDALGCP2_pred"` (spatial) or `"SDALGCP2_ST_pred"`
(spatio-temporal). For discrete spatial fits `$my_shp` is an `sf` with
`RR_mean`, `RR_se`, `ARR_mean`, `ARR_se` columns.
