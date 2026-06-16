# Predict relative risk from a fitted SDALGCP2 model

Predict relative risk from a fitted SDALGCP2 model

## Usage

``` r
# S3 method for class 'SDALGCP2'
predict(
  object,
  type = c("discrete", "continuous"),
  sampler = c("mcmc", "laplace"),
  cellsize = NULL,
  pred.loc = NULL,
  control.mcmc = NULL,
  ...
)
```

## Arguments

- object:

  an object of class `"SDALGCP2"` from
  [`SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2.md)
  or
  [`mcml_fit`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md).

- type:

  `"discrete"` for region-level inference or `"continuous"` for a
  spatially continuous surface.

- sampler:

  `"mcmc"` (MALA, default) or `"laplace"` (fast Gaussian approximation,
  no MCMC).

- cellsize:

  grid spacing for continuous prediction (ignored if `pred.loc`
  supplied).

- pred.loc:

  optional data frame of prediction coordinates (`x`, `y`) for
  continuous prediction.

- control.mcmc:

  optional MCMC controls; defaults to those used at fitting.

- ...:

  unused.

## Value

for `type = "discrete"`, an `sf` object augmented with posterior mean/SD
of incidence relative risk (`pMean_RR`/`pSD_RR`) and covariate-adjusted
relative risk (`pMean_ARR`/`pSD_ARR`); for `type = "continuous"`, a list
with the prediction grid and posterior summaries. Result carries class
`"SDALGCP2_pred"`.
