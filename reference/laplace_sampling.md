# Sample the latent field \[S \| Y\] (Poisson, non-nested) via C++ MALA

Fast drop-in for `SDALGCP::Laplace.sampling()` for the Poisson,
non-nested case: the Laplace mode (Newton) and the adaptive MALA loop
both run in C++. Given the same mode/covariance and seed it reproduces
the original R sampler bit-for-bit.

## Usage

``` r
laplace_sampling(mu, Sigma, y, units.m, control.mcmc)
```

## Arguments

- mu:

  prior mean vector.

- Sigma:

  prior covariance matrix.

- y:

  count vector.

- units.m:

  offset vector.

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md).

## Value

list with `samples` (kept x n matrix) and `h` (step sizes).
