# Map a fitted SDALGCP2 prediction

Maps a relative-risk surface (and its uncertainty) from
[`predict.SDALGCP2`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md).
Region-level predictions are drawn as choropleths; continuous
predictions as a raster, optionally masked to a boundary.

## Usage

``` r
# S3 method for class 'SDALGCP2_pred'
plot(x, variable = "RR", bound = NULL, midpoint = NULL, title = NULL, ...)
```

## Arguments

- x:

  an object of class `"SDALGCP2_pred"`.

- variable:

  for discrete: one of `"RR"`, `"SE_RR"`, `"ARR"`, `"SE_ARR"`; for
  continuous: `"RR"` or `"SE_RR"`.

- bound:

  optional `sf` boundary; continuous surfaces are masked to it and its
  outline overlaid.

- midpoint:

  optional value to centre a diverging colour scale (e.g. 1 for relative
  risk); if `NULL` a sequential viridis scale is used.

- title:

  optional plot title.

- ...:

  unused.

## Value

a `ggplot` object.
