# Map a fitted SDALGCP2 prediction

Maps any of the four predicted quantities from
[`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md)
– the relative risk `"relative_risk"`, the covariate-adjusted relative
risk `"adjusted_rr"`, or their standard errors
`"relative_risk_se"`/`"adjusted_rr_se"` – for either discrete
(choropleth) or continuous (raster) predictions.

## Usage

``` r
# S3 method for class 'SDALGCP2_pred'
plot(
  x,
  variable = c("relative_risk", "adjusted_rr", "relative_risk_se", "adjusted_rr_se"),
  bound = NULL,
  midpoint = NULL,
  title = NULL,
  ...
)
```

## Arguments

- x:

  an object of class `"SDALGCP2_pred"`.

- variable:

  one of `"relative_risk"`, `"adjusted_rr"`, `"relative_risk_se"`,
  `"adjusted_rr_se"`.

- bound:

  optional `sf` boundary; continuous surfaces are masked to it and its
  outline overlaid.

- midpoint:

  optional value to centre a diverging colour scale (defaults to 1 for
  the relative-risk columns, none for the standard errors).

- title:

  optional plot title.

- ...:

  unused.

## Value

a `ggplot` object.
