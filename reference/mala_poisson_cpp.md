# Adaptive MALA sampler for \[S \| Y\], Poisson, non-nested (C++)

Draw order per iteration is `d` normals then one uniform, giving
reproducible results under a common seed and the same mode/Sigma.tilde.

## Usage

``` r
mala_poisson_cpp(
  y,
  m,
  mu,
  Sigma,
  mode,
  Sigma_tilde,
  n_sim,
  burnin,
  thin,
  h_init,
  c1,
  c2
)
```

## Arguments

- y, m, mu, Sigma:

  data and prior as in `laplace_mode_poisson_cpp`.

- mode, Sigma_tilde:

  Laplace mode and covariance (preconditioner).

- n_sim, burnin, thin:

  MCMC length controls.

- h_init:

  initial step size; if not finite, `1.65 / d^(1/6)` is used.

- c1, c2:

  step-size adaptation constants.

## Value

list with `samples` (kept x d matrix of S draws) and `h`.
