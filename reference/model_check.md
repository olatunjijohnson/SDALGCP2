# Posterior-predictive model checking for an SDALGCP2 fit

Compares observed counts with fitted Poisson means, returns Pearson
residuals, and tests them for residual spatial autocorrelation with
Moran's I. A non-significant Moran's I indicates the spatial random
effect has absorbed the spatial structure.

## Usage

``` r
model_check(object, pred = NULL, nsim = 999, plot = TRUE)
```

## Arguments

- object:

  a fitted `"SDALGCP2"` object.

- pred:

  a discrete prediction from `predict(object, "discrete")`; if `NULL`
  one is computed with the fitting MCMC controls.

- nsim:

  permutations for the Moran's I p-value.

- plot:

  logical; draw the observed-vs-fitted scatter.

## Value

invisibly, a list with `fitted`, `residuals` and `moran`.
