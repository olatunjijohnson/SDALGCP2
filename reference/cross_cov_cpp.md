# Aggregated cross-covariance between prediction points and regions (C++)

Mirrors compute_cross_cov() for continuous prediction (bottleneck B5):
returns an n_pred x N matrix with entries sum_l w_jl \*
matern(\|\|x_pred - x_jl\|\|, phi, kappa) (weighted) or the mean.

## Usage

``` r
cross_cov_cpp(pred, coords, weights, phi, kappa, weighted, nthreads = 0L)
```

## Arguments

- pred:

  n_pred x 2 matrix of prediction coordinates.

- coords:

  list of N region point matrices.

- weights:

  list of N weight vectors, or empty for unweighted.

- phi:

  single spatial scale parameter.

- kappa:

  Matern smoothness.

- weighted:

  logical.

- nthreads:

  OpenMP threads (\<=0 default).
