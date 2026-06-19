# Sample the latent field \[S \| Y\] (Poisson, non-nested) via C++ MALA

Draws posterior samples of the latent Gaussian field for the Poisson,
non-nested case. The Laplace mode (Newton step) and the adaptive
Metropolis- adjusted Langevin (MALA) loop both run in C++ for speed,
with a fixed-seed path for reproducibility.

## Usage

``` r
laplace_sampling(mu, Sigma, y, units.m, control.mcmc)
```

## Arguments

- mu:

  prior mean vector.

- Sigma:

  prior covariance matrix.

- y:

  count vector.

- units.m:

  offset vector.

- control.mcmc:

  list from
  [`control_mcmc`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md).

## Value

list with `samples` (kept x n matrix) and `h` (step sizes).

## Examples

``` r
# \donttest{
## sample [S | Y] for a tiny 10-unit Poisson example
set.seed(1)
n <- 10
D <- as.matrix(dist(cbind(runif(n), runif(n))))
Sigma <- 0.4 * exp(-D / 0.3)
mu <- rep(log(2), n); m <- rep(100, n)
y <- rpois(n, m * exp(mu + as.numeric(t(chol(Sigma)) %*% rnorm(n))))
out <- laplace_sampling(mu, Sigma, y, m, control_mcmc(n.sim = 2000, burnin = 500, thin = 3))
dim(out$samples)     # (retained draws) x n
#> [1] 500  10
# }
```
