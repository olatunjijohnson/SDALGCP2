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

an object of class `"SDALGCP2_pred"`. For both `type`s it carries, for
every location (region or grid cell), the posterior mean and standard
error of two quantities:

- `RR` (relative risk):

  \\\exp(\eta)=\exp(d'\beta+S)\\ – the full relative risk, including the
  covariate effect.

- `ARR` (covariate-adjusted relative risk):

  \\\exp(S)\\ – the residual spatial relative risk after adjusting for
  covariates.

stored as `RR_mean`/`RR_se`/`ARR_mean`/`ARR_se` (and, for discrete fits,
as columns of the returned `sf`). Posterior draws are kept so that
[`exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md)
can be computed for either quantity.
