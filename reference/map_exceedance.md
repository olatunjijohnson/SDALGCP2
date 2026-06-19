# Map exceedance probabilities P(risk \> threshold)

Map exceedance probabilities P(risk \> threshold)

## Usage

``` r
map_exceedance(
  x,
  threshold = 1,
  which = c("adjusted_rr", "relative_risk"),
  bound = NULL,
  ...
)
```

## Arguments

- x:

  an `"SDALGCP2_pred"` object.

- threshold:

  a single relative-risk threshold.

- which:

  `"adjusted_rr"` (covariate-adjusted, default) or `"relative_risk"`.

- bound:

  optional `sf` boundary (continuous only).

- ...:

  unused.

## Value

a `ggplot` object.

## See also

[`exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md)
for the underlying probabilities.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
pr <- predict(fit, type = "discrete")
map_exceedance(pr, threshold = 1.5)           # P(adjusted RR > 1.5)

# }
```
