# Map a fitted SDALGCP2 prediction

Maps any of the four predicted quantities from
[`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md)
– relative risk `"RR"`, covariate-adjusted relative risk `"ARR"`, or
their standard errors `"RR_se"`/`"ARR_se"` – for either discrete
(choropleth) or continuous (raster) predictions.

## Usage

``` r
# S3 method for class 'SDALGCP2_pred'
plot(
  x,
  variable = c("RR", "ARR", "RR_se", "ARR_se"),
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

  one of `"RR"`, `"ARR"`, `"RR_se"`, `"ARR_se"`.

- bound:

  optional `sf` boundary; continuous surfaces are masked to it and its
  outline overlaid.

- midpoint:

  optional value to centre a diverging colour scale (defaults to 1 for
  `"RR"`/`"ARR"`, none for the standard errors).

- title:

  optional plot title.

- ...:

  unused.

## Value

a `ggplot` object.
