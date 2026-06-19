# Generate candidate sampling points inside each region

For every polygon feature in `my_shp` it produces candidate points and
aggregation weights, in the list format consumed by
[`precompute_corr`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md).

## Usage

``` r
sda_points(
  my_shp,
  delta,
  method = 1L,
  weighted = FALSE,
  pop_shp = NULL,
  rho = 0.55,
  giveup = 1000L
)
```

## Arguments

- my_shp:

  an `sf` object of `POLYGON`/`MULTIPOLYGON` features.

- delta:

  point spacing (grid step / SSI inhibition distance).

- method:

  1 = SSI (default), 2 = uniform random, 3 = regular grid.

- weighted:

  logical; if `TRUE`, weights are population density read from
  `pop_shp`, otherwise equal weights.

- pop_shp:

  a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of population density (required when `weighted = TRUE`).

- rho:

  packing density used to choose the number of points.

- giveup:

  SSI rejection limit.

## Value

a list of length `nrow(my_shp)`; each element has `xy` and `weight`.
Carries `"weighted"` and `"my_shp"` attributes.

## See also

[`precompute_corr`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md),
which consumes this output.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
pts <- sda_points(sdalgcp_data, delta = 1.2, method = 3)   # regular grid points
length(pts)          # one entry per region
#> [1] 64
str(pts[[1]])        # $xy candidate coordinates and $weight
#> List of 2
#>  $ xy    : num [1:9, 1:2] 0 1.2 2.4 0 1.2 2.4 0 1.2 2.4 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:9] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:2] "x" "y"
#>  $ weight: num [1:9] 0.111 0.111 0.111 0.111 0.111 ...
# }
```
