# Exceedance probabilities P(risk \> threshold)

Exceedance probabilities P(risk \> threshold)

## Usage

``` r
exceedance(object, thresholds, which = c("adjusted_rr", "relative_risk"))
```

## Arguments

- object:

  an `"SDALGCP2_pred"` object from
  [`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md).

- thresholds:

  numeric vector of thresholds.

- which:

  which quantity: `"adjusted_rr"` (the covariate-adjusted relative risk
  \\\exp(S)\\, default) or `"relative_risk"` (the relative risk
  \\\exp(d'\beta + S)\\).

## Value

a matrix of exceedance probabilities (locations x thresholds).

## See also

[`map_exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/map_exceedance.md)
to map them.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
pr <- predict(fit, type = "discrete")
## P(adjusted relative risk > 1) and > 1.5 for every region
ex <- exceedance(pr, thresholds = c(1, 1.5), which = "adjusted_rr")
head(ex)
#>           [,1]      [,2]
#> [1,] 0.4600000 0.1666667
#> [2,] 0.5400000 0.1600000
#> [3,] 0.3633333 0.0900000
#> [4,] 0.5633333 0.2533333
#> [5,] 0.9233333 0.6833333
#> [6,] 1.0000000 0.9866667
# }
```
