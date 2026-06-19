# Fit an SDA-LGCP with covariates measured on a different support

Covariates observed on a *different support* from the outcome (e.g. air-
quality monitors at point locations) are kriged to the candidate points
and enter the model on the intensity scale with a Berkson correction
that propagates the prediction uncertainty (see
`math/confounding-and-misalignment.pdf`).

## Usage

``` r
SDALGCP2_misaligned(
  formula,
  data,
  delta,
  covariates,
  phi = NULL,
  method = 3L,
  weighted = FALSE,
  pop_shp = NULL,
  berkson = TRUE,
  control.mcmc = NULL,
  max_iter = 10L,
  tol = 0.001,
  messages = FALSE
)
```

## Arguments

- formula:

  model formula; the covariate names appear on the right-hand side.

- data:

  `sf` polygons holding the response and offset (one row/region).

- delta:

  candidate-point spacing.

- covariates:

  a named list; each element is an `sf` carrying a column of the same
  name – the covariate's observed values on its own support, either
  **points** (e.g. monitors; plain kriging) or **polygons** (areal
  averages on a different partition; aggregated areal kriging).

- phi:

  spatial-scale grid for the outcome model (default from geometry).

- method, weighted, pop_shp:

  point-generation controls.

- berkson:

  logical; include the Berkson uncertainty correction (default `TRUE`).
  `FALSE` gives the naive kriged-mean plug-in.

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md).

- max_iter, tol:

  outer Gauss-Newton controls.

- messages:

  logical; print progress.

## Value

an object of class `"SDALGCP2"` with `misaligned = TRUE`.

## See also

[`SDALGCP2_raster`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_raster.md),
[`sdalgcp`](https://olatunjijohnson.github.io/SDALGCP2/reference/sdalgcp.md)

## Examples

``` r
# \donttest{
data(sdalgcp_data)
set.seed(1)
## a covariate z observed only at 40 scattered monitor points (a different support)
mon <- sf::st_as_sf(data.frame(x = runif(40, 0, 20), y = runif(40, 0, 20)),
                    coords = c("x", "y"))
mon$z <- scale(sf::st_coordinates(mon)[, 1])[, 1]
fit <- SDALGCP2_misaligned(cases ~ z + offset(log(pop)), sdalgcp_data, delta = 1.5,
                           covariates = list(z = mon),
                           control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
summary(fit)
#> Call: SDALGCP2_misaligned(formula = cases ~ z + offset(log(pop)), data = sdalgcp_data, 
#>     delta = 1.5, covariates = list(z = mon), control.mcmc = control_mcmc(n.sim = 2000, 
#>         burnin = 500, thin = 5))
#> 
#> Coefficients:
#>             Estimate Std.Err z value Pr(>|z|)    
#> (Intercept)   -6.154   0.194  -31.71  < 2e-16 ***
#> z              0.532   0.175    3.04   0.0023 ** 
#> sigma^2        0.805   0.190    4.25  2.2e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial scale phi: 2
#> Log-likelihood: 9.18732
#> MC importance-sampling ESS: 1 / 300 (0%);  log-lik MC SE: 0.925
#> Note: sigma^2 is the variance of the latent Gaussian process.
# }
```
