# Simulated aggregated disease-count data

A small, self-contained example dataset used throughout the help pages
and vignettes. It is simulated from the model the package fits: an 8x8
lattice of regions, a spatially structured covariate, a latent Gaussian
spatial field with exponential covariance, and Poisson counts with a
population offset. The true fixed effects are `(Intercept) = -6` and
`x1 = 0.6`; the latent field has variance \\\sigma^2 = 0.3\\ and
exponential scale \\\phi = 4\\.

## Usage

``` r
sdalgcp_data
```

## Format

An [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object of 64
`POLYGON` regions with columns:

- region:

  integer region identifier (1-64).

- cases:

  observed disease count in the region.

- x1:

  a standardised, spatially structured covariate.

- pop:

  population at risk (the offset; use `offset(log(pop))`).

- geometry:

  the region polygon.

## Source

Simulated; see `data-raw/sdalgcp_data.R` in the package sources.

## See also

[`liver`](https://olatunjijohnson.github.io/SDALGCP2/reference/liver.md)
for a real disease-count example.

## Examples

``` r
data(sdalgcp_data)
summary(sdalgcp_data$cases)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.000   1.000   4.000   6.031  10.000  25.000 
plot(sdalgcp_data["cases"])
```
