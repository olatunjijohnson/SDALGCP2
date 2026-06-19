# Monte Carlo maximum likelihood estimation for the spatial SDA-LGCP

Vectorised, Cholesky-based MCML estimation. Simulates the latent field
at an anchor, then profiles the importance-sampling MCML objective over
the supplied `phi` grid.

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

## See also

[`SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2.md)
(the end-to-end wrapper),
[`precompute_corr`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md)

## Examples

``` r
# \donttest{
data(sdalgcp_data)
df  <- sf::st_drop_geometry(sdalgcp_data)
pts <- sda_points(sdalgcp_data, delta = 1.2, method = 3)
cc  <- precompute_corr(pts, phi = seq(2, 8, length.out = 6))
fit <- mcml_fit(cases ~ x1 + offset(log(pop)), df, cc,
                control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
summary(fit)
#> Call: mcml_fit(formula = cases ~ x1 + offset(log(pop)), data = df, 
#>     corr = cc, control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, 
#>         thin = 5))
#> 
#> Coefficients:
#>             Estimate Std.Err z value Pr(>|z|)    
#> (Intercept)   -6.212   0.174  -35.67  < 2e-16 ***
#> x1             0.649   0.144    4.52  6.3e-06 ***
#> sigma^2        0.721   0.193    3.74  0.00018 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial scale phi: 2
#> Log-likelihood: 8.06238
#> MC importance-sampling ESS: 1 / 300 (0%);  log-lik MC SE: 0.953
#> Note: sigma^2 is the variance of the latent Gaussian process.
# }
```
