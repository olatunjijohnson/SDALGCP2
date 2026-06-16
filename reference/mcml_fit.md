# Monte Carlo maximum likelihood estimation for the spatial SDA-LGCP

Vectorised, Cholesky-based re-implementation of
`SDALGCP::SDALGCPParaEst()`. Simulates the latent field at an anchor,
then profiles the importance-sampling MCML objective over the supplied
`phi` grid.

## Usage

``` r
mcml_fit(
  formula,
  data,
  corr,
  par0 = NULL,
  control.mcmc = NULL,
  phi_method = c("grid", "direct"),
  nugget = FALSE,
  reanchor = 0L,
  reanchor_tol = 0.01,
  messages = FALSE
)
```

## Arguments

- formula:

  model formula, optionally with an
  [`offset()`](https://rdrr.io/r/stats/offset.html) term.

- data:

  data frame holding the model variables.

- corr:

  list with `R` (N x N x n_phi correlation array) and `phi`, e.g. from
  [`precompute_corr`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md).

- par0:

  optional starting values `c(beta, sigma2, phi)`; if `NULL` they are
  derived from a Poisson GLM.

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md)
  (defaults if `NULL`).

- phi_method:

  `"grid"` (profile over the `corr` phi grid, default) or `"direct"`
  (optimise phi continuously; exponential/Matern kernel).

- nugget:

  logical; if `TRUE` (requires `phi_method = "direct"`) add a relative
  nugget, fitting covariance \\\sigma^2(R(\phi)+\nu I)\\.

- reanchor:

  number of re-anchoring passes (re-simulate the latent field at the
  current optimum and refit) to raise the importance-sampling ESS.

- reanchor_tol:

  relative-change tolerance for stopping the re-anchoring loop.

- messages:

  logical; print optimiser progress.

## Value

an object of class `"SDALGCP2"` (estimates, covariance, profile, latent
samples and metadata).
