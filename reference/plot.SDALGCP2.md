# Plot an SDALGCP2 fit (the phi profile deviance)

Plot an SDALGCP2 fit (the phi profile deviance)

## Usage

``` r
# S3 method for class 'SDALGCP2'
plot(x, ...)
```

## Arguments

- x:

  an `"SDALGCP2"` object.

- ...:

  passed to
  [`phi_profile`](https://olatunjijohnson.github.io/SDALGCP2/reference/phi_profile.md).

## Value

invisibly, the profile (see
[`phi_profile`](https://olatunjijohnson.github.io/SDALGCP2/reference/phi_profile.md)).

## Examples

``` r
# \donttest{
data(sdalgcp_data)
fit <- SDALGCP2(cases ~ x1 + offset(log(pop)),
                sf::st_drop_geometry(sdalgcp_data), sdalgcp_data, delta = 1.2,
                control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
plot(fit)   # profile deviance for the spatial scale phi

# }
```
