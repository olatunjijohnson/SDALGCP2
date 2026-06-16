# Map an sdalgcp fit

Default visualisation: a choropleth of the covariate-adjusted relative
risk (spatial fits). Equivalent to `plot(predict(object), ...)`.

## Usage

``` r
# S3 method for class 'sdalgcp'
plot(
  x,
  type = c("risk", "incidence", "risk_se", "exceedance"),
  threshold = 1,
  ...
)
```

## Arguments

- x:

  an `"sdalgcp"` fit.

- type:

  `"risk"` (default), `"incidence"`, `"risk_se"` or `"exceedance"`.

- threshold:

  threshold for `type = "exceedance"`.

- ...:

  passed to the mapping layer.

## Value

a `ggplot` object.
