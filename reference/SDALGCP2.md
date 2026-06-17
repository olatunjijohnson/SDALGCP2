# Fit a spatial SDA-LGCP model

End-to-end user entry point: generates candidate points inside each
region, assembles the aggregated region-level correlation array (C++),
and estimates parameters by Monte Carlo maximum likelihood. This is the
modern, faster equivalent of `SDALGCP::SDALGCPMCML()` for the spatial
case.

## Usage

``` r
SDALGCP2(
  formula,
  data,
  my_shp,
  delta,
  phi = NULL,
  method = 1L,
  weighted = FALSE,
  pop_shp = NULL,
  kappa = 0.5,
  par0 = NULL,
  control.mcmc = NULL,
  phi_method = c("grid", "direct"),
  nugget = FALSE,
  confounding = c("none", "restricted"),
  reanchor = 0L,
  rho = 0.55,
  giveup = 1000L,
  nthreads = 0L,
  messages = FALSE
)
```

## Arguments

- formula:

  model formula, e.g. `cases ~ x1 + offset(log(pop))`.

- data:

  data frame with the model variables (one row per region).

- my_shp:

  `sf` polygons (or anything coercible via `st_as_sf`).

- delta:

  candidate-point spacing.

- phi:

  numeric vector of spatial scale parameters to profile; if `NULL`, a
  default grid from `sqrt(min area)` to `extent/10`.

- method:

  point method: 1 = SSI, 2 = uniform, 3 = regular grid.

- weighted:

  logical; population-weighted aggregation using `pop_shp`.

- pop_shp:

  population-density `SpatRaster` (needed if `weighted`).

- kappa:

  Matern smoothness for the spatial kernel (0.5 default).

- par0:

  optional starting values `c(beta, sigma2, phi)`.

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md).

- phi_method:

  how the spatial scale is estimated: `"grid"` (profile over the
  supplied `phi` grid, the robust default) or `"direct"` (optimise `phi`
  continuously inside the MCML objective; exponential kernel only). See
  the package vignette/PDF on the double-integral derivation.

- nugget:

  logical; if `TRUE` (requires `phi_method = "direct"`) add an
  unstructured region-level term, fitting covariance
  \\\sigma^2(R(\phi)+\nu I)\\ and estimating the relative nugget
  \\\nu=\tau^2/\sigma^2\\ with a standard error. Absorbs overdispersion
  not explained by the spatial structure.

- confounding:

  `"none"` (default) or `"restricted"` for restricted spatial regression
  (constrains the spatial random effect orthogonal to the fixed-effect
  design; fitted by a Laplace-approximate marginal likelihood).

- reanchor:

  number of re-anchoring passes: after fitting, the latent field is
  re-simulated at the current optimum and the model refit, which keeps
  the importance weights near-uniform (raises the MC effective sample
  size). 0 (default) fits once; 2-3 is usually ample.

- rho, giveup:

  point-generation controls.

- nthreads:

  OpenMP threads for the correlation build.

- messages:

  logical; print optimiser progress.

## Value

an object of class `"SDALGCP2"`.

## See also

[`mcml_fit`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md),
[`precompute_corr`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md),
[`sda_points`](https://olatunjijohnson.github.io/SDALGCP2/reference/sda_points.md)

## Examples

``` r
# \donttest{
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
## ---- simulate a lattice of regions and aggregated counts ----
set.seed(1)
bound <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20)))
shp   <- st_sf(geometry = st_make_grid(bound, n = c(8, 8)))
N     <- nrow(shp)

pts   <- sda_points(shp, delta = 1.2, method = 3)      # regular grid points
phi_grid <- seq(1, 5, length.out = 8)
corr  <- precompute_corr(pts, phi_grid)
Sig   <- 0.5 * corr$R[, , which.min(abs(phi_grid - 2.5))]
x1    <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
pop   <- round(runif(N, 500, 3000))
y     <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
                            as.numeric(t(chol(Sig)) %*% rnorm(N))))
dat   <- data.frame(y = y, x1 = x1, pop = pop)

## ---- fit ----
ctrl <- control_mcmc(n.sim = 6000, burnin = 1500, thin = 6, h = 1.65 / N^(1/6))
fit  <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.2,
                 phi = phi_grid, method = 3, control.mcmc = ctrl)
summary(fit)
#> Call: SDALGCP2(formula = y ~ x1 + offset(log(pop)), data = dat, my_shp = shp, 
#>     delta = 1.2, phi = phi_grid, method = 3, control.mcmc = ctrl)
#> 
#> Coefficients:
#>             Estimate Std.Err z value Pr(>|z|)    
#> (Intercept)   -5.966   0.125  -47.62  < 2e-16 ***
#> x1             0.824   0.107    7.72  1.2e-14 ***
#> sigma^2        0.510   0.207    2.47    0.014 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial scale phi: 1.57143
#> Log-likelihood: 8.33903
#> MC importance-sampling ESS: 1 / 750 (0%);  log-lik MC SE: 0.876
#> Note: sigma^2 is the variance of the latent Gaussian process.

## ---- predict ----
pred_d <- predict(fit, type = "discrete",   sampler = "mcmc",    control.mcmc = ctrl)
pred_c <- predict(fit, type = "continuous", sampler = "laplace", cellsize = 1,
                  control.mcmc = ctrl)
#> Error in predict.SDALGCP2(fit, type = "continuous", sampler = "laplace",     cellsize = 1, control.mcmc = ctrl): 'list' object cannot be coerced to type 'integer'
# }
```
