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
