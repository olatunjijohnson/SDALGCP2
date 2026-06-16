# Aggregated correlation and its phi-derivatives at one phi (C++, Matern)

Aggregated correlation and its phi-derivatives at one phi (C++, Matern)

## Usage

``` r
corr_and_grad_cpp(coords, weights, phi, kappa, weighted, nthreads = 0L)
```

## Arguments

- coords:

  list of N candidate-point matrices (n_i x 2).

- weights:

  list of N weight vectors (each summing to 1), or empty for the
  unweighted (mean) case.

- phi:

  spatial scale (\> 0).

- kappa:

  Matern smoothness; one of 0.5, 1.5, 2.5.

- weighted:

  logical; population-weighted aggregation.

- nthreads:

  OpenMP threads (\<= 0 = default).

## Value

list with N x N matrices `R`, `dR` (dR/dphi) and `d2R` (d2R/dphi2).
