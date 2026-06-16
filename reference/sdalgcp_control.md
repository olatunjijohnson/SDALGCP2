# Control settings for [`sdalgcp`](https://olatunjijohnson.github.io/SDALGCP2/reference/sdalgcp.md)

Bundles the technical knobs so that a default fit needs none of them.

## Usage

``` r
sdalgcp_control(
  delta = NULL,
  points_per_region = 16,
  point_method = c("regular", "uniform", "ssi"),
  scale = c("continuous", "grid"),
  phi = NULL,
  kappa = 0.5,
  kappa_t = 0.5,
  nugget = FALSE,
  confounding = c("none", "restricted"),
  reanchor = 2L,
  n_sim = 10000L,
  burnin = 2000L,
  thin = 8L,
  tilt_spatial = FALSE,
  nthreads = 0L
)
```

## Arguments

- delta:

  candidate-point spacing. If `NULL` (default) it is chosen
  automatically to place roughly `points_per_region` points in a typical
  region.

- points_per_region:

  target number of candidate points per region used to pick `delta`
  automatically.

- point_method:

  how candidate points are laid out: `"regular"` (deterministic grid,
  default), `"uniform"` or `"ssi"`.

- scale:

  how the spatial scale \\\phi\\ is estimated: `"continuous"` (optimised
  directly, no grid – the default) or `"grid"` (profiled over `phi`).
  Spatio-temporal fits always profile \\\phi\\ on a grid.

- phi:

  optional \\\phi\\ grid (only used when `scale = "grid"` or for
  spatio-temporal fits); chosen from the geometry if `NULL`.

- kappa:

  spatial Matern smoothness (`0.5`, `1.5` or `2.5`).

- kappa_t:

  temporal Matern smoothness (spatio-temporal fits).

- nugget:

  logical; add an unstructured region-level term (overdispersion).
  Requires `scale = "continuous"`.

- confounding:

  `"none"` (default) or `"restricted"`. With `"restricted"`, restricted
  spatial regression is used: the spatial random effect is constrained
  to the orthogonal complement of the fixed-effect design so it cannot
  absorb a spatially structured covariate (avoids spatial confounding /
  attenuation of `beta`). Spatial models only.

- reanchor:

  number of re-anchoring passes (re-simulate the latent field at the
  optimum and refit) for reliable variance estimates. Default `2`.

- n_sim, burnin, thin:

  MCMC length controls for the latent-field sampler.

- tilt_spatial:

  logical; for raster covariates, use the fully covariate-tilted
  correlation (see
  [`SDALGCP2_raster`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_raster.md)).

- nthreads:

  OpenMP threads for the correlation assembly (0 = default).

## Value

a list of control settings.
