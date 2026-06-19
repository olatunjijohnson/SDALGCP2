# Coefficient plot of fixed effects (and sigma^2) with confidence intervals

Coefficient plot of fixed effects (and sigma^2) with confidence
intervals

## Usage

``` r
coef_plot(object, level = 0.95, intercept = FALSE)
```

## Arguments

- object:

  a fitted `"SDALGCP2"` object.

- level:

  confidence level.

- intercept:

  logical; include the intercept.

## Value

a `ggplot` object.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
coef_plot(fit)
#> `height` was translated to `width`.

# }
```
