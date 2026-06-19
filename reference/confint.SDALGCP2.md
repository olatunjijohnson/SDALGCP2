# Wald confidence intervals for an SDALGCP2 fit

Wald confidence intervals for an SDALGCP2 fit

## Usage

``` r
# S3 method for class 'SDALGCP2'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  an object of class `"SDALGCP2"`.

- parm:

  parameters to report (names or indices); default all.

- level:

  confidence level.

- ...:

  unused.

## Value

a matrix of lower/upper confidence limits.

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
               control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
                                         reanchor = 0))
confint(fit)
#>                   2.5%     97.5%
#> (Intercept) -6.5712603 -5.942619
#> x1           0.4248018  0.973310
#> sigma^2      0.6204925  1.674251
#> phi          0.4610345  2.003098
# }
```
