# Predict relative risk from an sdalgcp fit

Returns a prediction object carrying, for every location, the posterior
mean and standard error of the relative risk `relative_risk`
(\\\exp(\eta)=\exp(d'\beta+S)\\) and the covariate-adjusted relative
risk `adjusted_rr` (\\\exp(S)\\). Map it with
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

for a spatial fit, an `sf` of class `"SDALGCP2_pred"` with
`relative_risk`, `relative_risk_se`, `adjusted_rr` and `adjusted_rr_se`
columns (polygons for `type = "discrete"`, grid points for
`"continuous"`); for a spatio-temporal fit, an `"SDALGCP2_ST_pred"`
object (see
[`predict.SDALGCP2_ST`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2_ST.md)).
