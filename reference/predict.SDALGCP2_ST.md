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

a list with \\N\times T\\ matrices `RR_mean`, `RR_sd`, `ARR_mean`,
`ARR_sd` and a long data frame `table`.
