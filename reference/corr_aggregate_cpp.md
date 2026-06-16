# Aggregated correlation array (C++)

Aggregated correlation array (C++)

## Usage

``` r
corr_aggregate_cpp(coords, weights, phi, kappa, weighted, nthreads = 0L)
```

## Arguments

- coords:

  list of length N; element i is an n_i x 2 numeric matrix of
  candidate-point coordinates inside region i.

- weights:

  list of length N of weight vectors (each summing to 1), or an empty
  list for the unweighted (mean) case.

- phi:

  numeric vector of spatial scale parameters.

- kappa:

  Matern smoothness (0.5, 1.5 or 2.5 use closed forms).

- weighted:

  logical; TRUE for population-weighted aggregation.

- nthreads:

  number of OpenMP threads (\<=0 uses the OpenMP default).

## Value

a numeric array of dimension N x N x length(phi).
