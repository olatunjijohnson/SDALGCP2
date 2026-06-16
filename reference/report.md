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
