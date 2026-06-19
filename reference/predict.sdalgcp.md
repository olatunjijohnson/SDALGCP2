# Predict relative risk from an sdalgcp fit

Returns a prediction object carrying, for every location, the posterior
mean and standard error of the relative risk `relative_risk`
(\\\exp(\eta)=\exp(d'\beta+S)\\) and the covariate-adjusted relative
risk `adjusted_rr` (\\\exp(S)\\). Map it with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and get hotspot
probabilities with
[`exceedance`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md).

## Usage

``` r
# S3 method for class 'sdalgcp'
predict(
  object,
  type = c("discrete", "continuous"),
  sampler = c("mcmc", "laplace"),
  cellsize = NULL,
  ...
)
```

## Arguments

- object:

  an `"sdalgcp"` fit.

- type:

  `"discrete"` (region level, default) or `"continuous"` (a grid
  surface). Ignored for spatio-temporal fits.

- sampler:

  `"mcmc"` (default) or `"laplace"`.

- cellsize:

  grid spacing for `type = "continuous"`.

- ...:

  passed to the underlying predictor.

## Value

for a spatial fit, an `sf` of class `"SDALGCP2_pred"` with
`relative_risk`, `relative_risk_se`, `adjusted_rr` and `adjusted_rr_se`
columns (polygons for `type = "discrete"`, grid points for
`"continuous"`); for a spatio-temporal fit, an `"SDALGCP2_ST_pred"`
object (see
[`predict.SDALGCP2_ST`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2_ST.md)).

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
pr <- predict(fit)            # discrete by default; an sf of relative risks
head(pr)
#> Simple feature collection with 6 features and 8 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 0 xmax: 15 ymax: 2.5
#> CRS:           NA
#>   region cases          x1  pop                       geometry relative_risk
#> 1      1     2 -2.03331626 3840 POLYGON ((0 0, 2.5 0, 2.5 2...  0.0005317242
#> 2      2     3 -1.64601792 3985 POLYGON ((2.5 0, 5 0, 5 2.5...  0.0007289839
#> 3      3     1 -1.25871959 2236 POLYGON ((5 0, 7.5 0, 7.5 2...  0.0008197207
#> 4      4     0 -0.87142125  846 POLYGON ((7.5 0, 10 0, 10 2...  0.0012181931
#> 5      5     3 -0.48412292  874 POLYGON ((10 0, 12.5 0, 12....  0.0025647230
#> 6      6    12 -0.09682458 2231 POLYGON ((12.5 0, 15 0, 15 ...  0.0046964955
#>   relative_risk_se adjusted_rr adjusted_rr_se
#> 1     0.0002366949   0.9517023      0.4236464
#> 2     0.0003048683   1.0135720      0.4238859
#> 3     0.0003512447   0.8853696      0.3793747
#> 4     0.0005202427   1.0221081      0.4365024
#> 5     0.0010285486   1.6716415      0.6703899
#> 6     0.0012838339   2.3779273      0.6500302
# }
```
