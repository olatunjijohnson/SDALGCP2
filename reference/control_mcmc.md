# MCMC control settings for the MALA sampler

MCMC control settings for the MALA sampler

## Usage

``` r
control_mcmc(
  n.sim = 10000,
  burnin = 2000,
  thin = 8,
  h = NULL,
  c1.h = 0.01,
  c2.h = 1e-04
)
```

## Arguments

- n.sim:

  total number of iterations.

- burnin:

  burn-in iterations to discard.

- thin:

  thinning interval; `(n.sim - burnin)` must be a multiple.

- h:

  initial Langevin step size; if missing, `1.65 / d^(1/6)` is used.

- c1.h, c2.h:

  step-size adaptation constants.

## Value

a named list consumed by
[`laplace_sampling`](https://olatunjijohnson.github.io/SDALGCP2/reference/laplace_sampling.md)
/ the fit.
