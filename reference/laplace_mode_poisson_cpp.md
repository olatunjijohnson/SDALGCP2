# Conditional mode and Laplace covariance for \[S \| Y\], Poisson, non-nested (C++)

Conditional mode and Laplace covariance for \[S \| Y\], Poisson,
non-nested (C++)

## Usage

``` r
laplace_mode_poisson_cpp(y, m, mu, Sigma, tol = 1e-08, maxit = 100L)
```

## Arguments

- y:

  count vector.

- m:

  offset vector (e.g. expected counts / population).

- mu:

  prior mean vector of the latent field.

- Sigma:

  prior covariance matrix.

- tol:

  convergence tolerance on the gradient infinity-norm.

- maxit:

  maximum Newton iterations.

## Value

list with `mode` and `Sigma_tilde`.
