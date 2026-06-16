# Reference (pure-R) aggregated correlation builder

Slow but dependency-light implementation kept for correctness testing
and benchmarking against
[`precompute_corr`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md).
Computes the exponential (`kappa = 0.5`) aggregated correlation only.

## Usage

``` r
precompute_corr_ref(points, phi, weighted = NULL)
```

## Arguments

- points:

  a list of length \\N\\; each element holds `$xy` (an \\n_i \times 2\\
  matrix of candidate-point coordinates) and, when weighted, `$weight`
  (a length-\\n_i\\ vector summing to 1). The `"weighted"` and
  `"my_shp"` attributes produced by the point-generation step are
  honoured and carried through.

- phi:

  numeric vector of spatial scale parameters.

- weighted:

  logical; if `NULL` (default) it is taken from
  `attr(points, "weighted")`.

## Value

the \\N \times N \times\\ `length(phi)` correlation array.
