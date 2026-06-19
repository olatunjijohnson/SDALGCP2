# Profile likelihood and confidence interval for the spatial scale phi

Spline-smoothed profile deviance for `phi`, with the `coverage`-level
confidence interval where the deviance crosses the chi-squared cutoff.

## Usage

``` r
phi_profile(object, coverage = 0.95, plot = TRUE)
```

## Arguments

- object:

  a fitted `"SDALGCP2"` object.

- coverage:

  confidence level.

- plot:

  logical; draw the deviance curve.

## Value

invisibly, a list with the interval and the smoothed profile; a `ggplot`
is drawn when `plot = TRUE`.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
## profile phi on a grid (scale = "grid") so there is a deviance curve to draw
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(scale = "grid", n_sim = 2000,
                                         burnin = 500, thin = 5, reanchor = 0))
phi_profile(fit)

# }
```
