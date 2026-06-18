# Discrete (region x time) prediction for a spatio-temporal fit

Draws the latent field at the fitted optimum and returns posterior mean
and SD of the incidence relative risk \\\exp(\mu+S)\\ and
covariate-adjusted relative risk \\\exp(S)\\ for every region and time.

## Usage

``` r
# S3 method for class 'SDALGCP2_ST'
predict(object, control.mcmc = NULL, ...)
```

## Arguments

- object:

  an `"SDALGCP2_ST"` fit.

- control.mcmc:

  optional MCMC controls (defaults to the fitting ones).

- ...:

  unused.

## Value

an `"SDALGCP2_ST_pred"` object with \\N\times T\\ matrices
`relative_risk`, `relative_risk_se` (\\\exp(\mu+S)\\) and `adjusted_rr`,
`adjusted_rr_se` (\\\exp(S)\\), a long `table`, and the geometry; map a
time slice with
[`plot.SDALGCP2_ST_pred`](https://olatunjijohnson.github.io/SDALGCP2/reference/plot.SDALGCP2_ST_pred.md).
