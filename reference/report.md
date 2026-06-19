# One-call panel of post-fit graphics

Returns the maps and summaries an analyst usually wants after fitting:
relative-risk and uncertainty maps, an exceedance map, the coefficient
plot and the phi profile. The pieces are returned as a named list of
`ggplot` objects so they can be arranged or printed individually.

## Usage

``` r
report(object, pred = NULL, threshold = 1.5, ...)
```

## Arguments

- object:

  a fitted `"SDALGCP2"` object.

- pred:

  optional discrete prediction; computed if `NULL`.

- threshold:

  relative-risk threshold for the exceedance map.

- ...:

  passed to
  [`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md)
  when `pred` is computed.

## Value

a named list of `ggplot` objects.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
figs <- report(fit, threshold = 1.5)
names(figs)          # relative_risk, uncertainty, exceedance, coefficients, ...
#> [1] "relative_risk" "uncertainty"   "exceedance"    "coefficients" 
#> [5] "phi_profile"  
figs$relative_risk   # print one of the maps

# }
```
