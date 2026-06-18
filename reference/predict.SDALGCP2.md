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

an `sf` (class `c("SDALGCP2_pred", "sf", "data.frame")`) with one row
per location – polygons for `type = "discrete"`, grid-cell points for
`type = "continuous"` – carrying the posterior mean and standard error
of two relative-risk quantities:

- `relative_risk`, `relative_risk_se`:

  the relative risk \\\exp(d'\beta + S)\\ – the fitted risk relative to
  the offset baseline, combining the covariate effect and the residual
  spatial variation. This is the headline disease-mapping quantity.

- `adjusted_rr`, `adjusted_rr_se`:

  the covariate-adjusted relative risk \\\exp(S)\\ – the purely spatial
  relative risk that remains after holding the covariates fixed (the
  spatial signal the covariates do not explain).

The full posterior draws are retained as object attributes so that
[`exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md)
and
[`map_exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/map_exceedance.md)
can be computed for either quantity. Map a column with
[`plot.SDALGCP2_pred`](https://olatunjijohnson.github.io/SDALGCP2/reference/plot.SDALGCP2_pred.md).
