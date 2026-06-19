# Fit a spatio-temporal SDA-LGCP model (Kronecker-free)

Separable space-time SDA-LGCP for aggregated counts observed over the
same `N` regions at `T` times. The spatial scale `phi` is profiled on a
grid; the temporal Matern range `nu` is estimated continuously. The
likelihood never forms the \\(NT)\times(NT)\\ covariance.

## Usage

``` r
SDALGCP2_ST(
  formula,
  data,
  my_shp,
  times,
  delta,
  phi = NULL,
  kappa = 0.5,
  kappa_t = 0.5,
  method = 3L,
  weighted = FALSE,
  pop_shp = NULL,
  control.mcmc = NULL,
  reanchor = 0L,
  messages = FALSE
)
```

## Arguments

- formula:

  model formula (with optional `offset(log(pop))`).

- data:

  data frame of `N*T` rows ordered by time then region (rows
  `(t-1)*N + 1:N` are time `t`).

- my_shp:

  `sf` polygons for the `N` regions.

- times:

  numeric vector of length `T` of observation times.

- delta:

  candidate-point spacing.

- phi:

  spatial-scale grid (default from geometry).

- kappa:

  spatial Matern smoothness.

- kappa_t:

  temporal Matern smoothness.

- method, weighted, pop_shp:

  point-generation controls.

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md).

- reanchor:

  number of re-anchoring passes (re-simulate the latent field at the
  current optimum and refit); improves the variance-parameter estimates.

- messages:

  logical; print progress.

## Value

an object of class `c("SDALGCP2_ST","SDALGCP2")` with `phi_opt`,
`nu_opt`, coefficient table and covariance.

## See also

[`sdalgcp`](https://olatunjijohnson.github.io/SDALGCP2/reference/sdalgcp.md)
(friendly wrapper),
[`predict.SDALGCP2_ST`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2_ST.md)

## Examples

``` r
# \donttest{
data(sdalgcp_data)
shp <- sdalgcp_data
## build a 3-time panel (data frame, N*T rows ordered by time then region)
times <- 1:3
dat <- do.call(rbind, lapply(times, function(t) {
  d <- sf::st_drop_geometry(shp); d$time <- t
  d$cases <- rpois(nrow(d), d$pop * exp(-6 + 0.6 * d$x1 + 0.1 * (t - 2)))
  d
}))
fit <- SDALGCP2_ST(cases ~ x1 + offset(log(pop)), dat, shp, times = times,
                   delta = 1.5,
                   control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
fit$phi_opt; fit$nu_opt
#> [1] 2
#>           
#> 0.8773136 
# }
```
