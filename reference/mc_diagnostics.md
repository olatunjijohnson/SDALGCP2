# Importance-sampling diagnostics for an MCML fit

The MCML estimate reweights latent samples drawn at the anchor towards
the optimum. When the optimum is far from the anchor the weights become
uneven and the estimate unreliable. This reports the effective sample
size of the importance weights at the maximiser and a Monte Carlo
standard error for the maximised log-likelihood,
\\\mathrm{SE}\approx\sqrt{1/\mathrm{ESS}-1/B}\\.

## Usage

``` r
mc_diagnostics(object, warn_frac = 0.1)
```

## Arguments

- object:

  a fitted `"SDALGCP2"` object.

- warn_frac:

  warn if the ESS falls below this fraction of \\B\\.

## Value

invisibly, a list with `B`, `ESS`, `ESS_frac` and `se_loglik`.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
d <- mc_diagnostics(fit)
#> Warning: Low importance-sampling ESS (1 of 300, 0.4%): consider re-anchoring (iterate = TRUE) or a par0 closer to the optimum.
d$ESS_frac           # importance-sampling ESS as a fraction of the draws
#> [1] 0.003682337
# }
```
