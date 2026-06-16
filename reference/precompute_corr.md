# Precompute aggregated region-level correlation matrices

Builds the \\N \times N \times\\ `length(phi)` array of region-level
correlations used by the SDA-LGCP model, where \$\$R(\phi)\_{ij} =
\sum\_{k,l} w\_{ik} w\_{jl}\\ C(\lVert x\_{ik}-x\_{jl}\rVert; \phi,
\kappa)\$\$ (population-weighted) or the unweighted mean over
candidate-point pairs. The heavy reduction runs in C++ (OpenMP-parallel
over region pairs); this is the modern, fast replacement for
`SDALGCP::precomputeCorrMatrix()`.

## Usage

``` r
precompute_corr(points, phi, kappa = 0.5, weighted = NULL, nthreads = 0L)
```

## Arguments

- points:

  a list of length \\N\\; each element holds `$xy` (an \\n_i \times 2\\
  matrix of candidate-point coordinates) and, when weighted, `$weight`
  (a length-\\n_i\\ vector summing to 1). The `"weighted"` and
  `"my_shp"` attributes produced by the point-generation step are
  honoured and carried through.

- phi:

  numeric vector of spatial scale parameters.

- kappa:

  Matern smoothness; `0.5` (exponential, default), `1.5` or `2.5` use
  closed forms in C++.

- weighted:

  logical; if `NULL` (default) it is taken from
  `attr(points, "weighted")`.

- nthreads:

  number of OpenMP threads; `0` (default) uses the OpenMP runtime
  default.

## Value

a list with `R` (the correlation array) and `phi`, carrying `weighted`,
`my_shp` and `S_coord` attributes on `R`.
