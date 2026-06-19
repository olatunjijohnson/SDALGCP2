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
  rasters = NULL,
  covariates = NULL,
  confounding = c("none", "restricted"),
  berkson = TRUE,
  max_iter = 10L,
  tol = 0.001,
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

- rasters:

  optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of spatially continuous, time-invariant covariates (layers named in
  `formula`); they enter on the intensity scale as in
  [`SDALGCP2_raster`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_raster.md),
  fitted by a Gauss-Newton tilting loop around the space-time
  likelihood.

- covariates:

  optional named list of `sf` covariate layers measured on a different
  (time-invariant) support; each is kriged to the candidate points with
  a Berkson correction as in
  [`SDALGCP2_misaligned`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_misaligned.md).

- confounding:

  `"none"` (default) or `"restricted"` for restricted spatial regression
  against space-time confounding (see Details).

- berkson:

  logical; include the Berkson uncertainty correction for `covariates`
  (default `TRUE`).

- max_iter, tol:

  Gauss-Newton controls for the `rasters`/`covariates` tilting loop.

- messages:

  logical; print progress.

## Value

an object of class `c("SDALGCP2_ST","SDALGCP2")` with `phi_opt`,
`nu_opt`, coefficient table and covariance.

## Details

With `rasters` or `covariates` the covariate surface is taken to be
constant over time (time-varying covariates can still be supplied as
ordinary columns of `data`). `confounding = "restricted"` constrains the
space-time random effect to the orthogonal complement of the
fixed-effect design and is fitted by an analytic Laplace-marginal
likelihood; it reduces to the spatial restricted fit when `T = 1` and is
not currently combined with `rasters`/`covariates`.

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

## restricted spatial regression against space-time confounding
fit_c <- SDALGCP2_ST(cases ~ x1 + offset(log(pop)), dat, shp, times = times,
                     delta = 1.5, phi = c(2, 4, 6), confounding = "restricted")
fit_c$beta_opt
#> (Intercept)          x1 
#>  -5.9945200   0.5845133 

## a spatially continuous (raster) covariate, aggregated on the intensity scale
r <- terra::rast(terra::ext(0, 20, 0, 20), resolution = 0.5)
terra::values(r) <- as.numeric(scale(terra::crds(r)[, 1])); names(r) <- "z"
fit_r <- SDALGCP2_ST(cases ~ z + offset(log(pop)), dat, shp, times = times,
                     delta = 1.5, phi = c(2, 4, 6), rasters = r, max_iter = 4,
                     control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
fit_r$beta_opt
#> (Intercept)           z 
#>   -6.005052    0.490473 
# }
```
