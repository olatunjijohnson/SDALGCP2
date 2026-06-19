# Fit an SDA-LGCP with spatially continuous (raster) covariates

Covariates supplied as rasters enter the model at the candidate-point
level and are aggregated on the intensity (exp) scale via a log-sum-exp
offset \\b_i(\beta)=\log\sum_k w\_{ik}\exp(z(x\_{ik})^\top\beta)\\ – the
statistically correct alternative to averaging the predictor over each
polygon. Estimation is a Gauss-Newton fixed point that reuses
[`mcml_fit`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md)
with the intensity-tilted effective design.

## Usage

``` r
SDALGCP2_raster(
  formula,
  data,
  my_shp,
  delta,
  rasters,
  phi = NULL,
  method = 3L,
  weighted = FALSE,
  pop_shp = NULL,
  kappa = 0.5,
  tilt_spatial = FALSE,
  control.mcmc = NULL,
  max_iter = 10L,
  tol = 0.001,
  messages = FALSE
)
```

## Arguments

- formula:

  model formula; right-hand-side names must match raster layer names.
  The response and an `offset(log(pop))` come from `data`.

- data:

  data frame with the response and offset (one row per region).

- my_shp:

  `sf` polygons.

- delta:

  candidate-point spacing.

- rasters:

  a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  (or object coercible by
  [`terra::rast`](https://rspatial.github.io/terra/reference/rast.html))
  whose layers are the spatially varying covariates named in `formula`.

- phi:

  spatial-scale grid (default chosen from the geometry).

- method, weighted, pop_shp:

  point-generation controls (see
  [`sda_points`](https://olatunjijohnson.github.io/SDALGCP2/reference/sda_points.md)).

- kappa:

  Matern smoothness for the spatial correlation.

- tilt_spatial:

  logical; if `FALSE` (default) the spatial correlation uses the
  population weights and is precomputed once (covariates enter only
  through the log-sum-exp offset). If `TRUE`, the correlation
  \\R^c(\beta)\\ is rebuilt each iteration from the intensity-tilted
  weights \\c\_{ik}(\beta)\\ and a log-normal aggregation correction
  \\\tfrac12\sigma^2(1-R^c\_{ii})\\ is added – the fully tilted model
  (more accurate, more costly).

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md).

- max_iter, tol:

  outer Gauss-Newton iteration controls.

- messages:

  logical; print progress.

## Value

an object of class `"SDALGCP2"` (as
[`mcml_fit`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md))
with extra fields `raster = TRUE` and `n_iter`.

## See also

[`SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2.md),
[`mcml_fit`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md)

## Examples

``` r
# \donttest{
data(sdalgcp_data)
## a spatially continuous covariate supplied as a raster layer named "z"
r <- terra::rast(terra::ext(0, 20, 0, 20), resolution = 0.5)
terra::values(r) <- as.numeric(scale(terra::crds(r)[, 1]))   # west-east gradient
names(r) <- "z"
df <- sf::st_drop_geometry(sdalgcp_data)
fit <- SDALGCP2_raster(cases ~ z + offset(log(pop)), df, sdalgcp_data,
                       delta = 1.5, rasters = r,
                       control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
summary(fit)
#> Call: SDALGCP2_raster(formula = cases ~ z + offset(log(pop)), data = df, 
#>     my_shp = sdalgcp_data, delta = 1.5, rasters = r, control.mcmc = control_mcmc(n.sim = 2000, 
#>         burnin = 500, thin = 5))
#> 
#> Coefficients:
#>             Estimate Std.Err z value Pr(>|z|)    
#> (Intercept)   -6.187   0.210  -29.51  < 2e-16 ***
#> z              0.676   0.171    3.95  7.9e-05 ***
#> sigma^2        0.891   0.209    4.27  2.0e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial scale phi: 2
#> Log-likelihood: 7.8569
#> MC importance-sampling ESS: 3 / 300 (1%);  log-lik MC SE: 0.63
#> Note: sigma^2 is the variance of the latent Gaussian process.
# }
```
