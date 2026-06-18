# Generates data/sdalgcp_data.rda: a small, self-contained simulated disease-count
# dataset used throughout the examples and tutorials. Re-run with:
#   Rscript data-raw/sdalgcp_data.R
#
# The data are simulated from the same generative model the package fits: an
# 8x8 lattice of regions, a spatially structured covariate x1, a latent Gaussian
# field S with exponential covariance over the region centroids, and Poisson
# counts with a population offset and true coefficients (Intercept = -6, x1 = 0.6).

suppressMessages(library(sf))

set.seed(2024)

bound   <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20)))
grid    <- st_make_grid(bound, n = c(8, 8))
regions <- st_sf(region = seq_along(grid), geometry = grid)

ctr <- st_coordinates(st_centroid(regions))           # region centroids

## spatially structured covariate (a SW-NE gradient), standardised
x1 <- as.numeric(scale(ctr[, 1] + 0.5 * ctr[, 2]))

## latent spatial field: exponential covariance over centroids
D      <- as.matrix(dist(ctr))
phi    <- 4
sigma2 <- 0.3
Sigma  <- sigma2 * exp(-D / phi)
S      <- as.numeric(t(chol(Sigma)) %*% rnorm(nrow(regions)))

## population offset and Poisson counts
pop   <- round(runif(nrow(regions), 800, 4000))
eta   <- -6 + 0.6 * x1 + S                              # true beta = c(-6, 0.6)
cases <- rpois(nrow(regions), pop * exp(eta))

sdalgcp_data <- regions
sdalgcp_data$x1    <- x1
sdalgcp_data$pop   <- pop
sdalgcp_data$cases <- cases
sdalgcp_data <- sdalgcp_data[, c("region", "cases", "x1", "pop")]

save(sdalgcp_data, file = "data/sdalgcp_data.rda", compress = "xz")
