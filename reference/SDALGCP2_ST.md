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
