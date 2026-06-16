# Predict relative risk from an sdalgcp fit

Returns the fitted region-level relative risk as an `sf` object (for
spatial fits) so it can be mapped directly, or a long data frame for
spatio-temporal fits.

## Usage

``` r
# S3 method for class 'sdalgcp'
predict(
  object,
  type = c("risk", "incidence", "exceedance"),
  threshold = 1,
  ...
)
```

## Arguments

- object:

  an `"sdalgcp"` fit.

- type:

  `"risk"` for covariate-adjusted relative risk \\\exp(S)\\ (default),
  `"incidence"` for \\\exp(\mu+S)\\, or `"exceedance"` for
  \\P(\mathrm{risk} \> \mathrm{threshold})\\.

- threshold:

  threshold for `type = "exceedance"`.

- ...:

  passed to the underlying predictor.

## Value

for spatial fits, the model's `sf` augmented with `relative_risk`,
`relative_risk_se` (and `incidence`, `exceedance` as requested); for
spatio-temporal fits, a list with region-by-time matrices and a long
table.
