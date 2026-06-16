#' Fit a spatial SDA-LGCP model
#'
#' End-to-end user entry point: generates candidate points inside each region,
#' assembles the aggregated region-level correlation array (C++), and estimates
#' parameters by Monte Carlo maximum likelihood. This is the modern, faster
#' equivalent of \code{SDALGCP::SDALGCPMCML()} for the spatial case.
#'
#' @param formula model formula, e.g. \code{cases ~ x1 + offset(log(pop))}.
#' @param data data frame with the model variables (one row per region).
#' @param my_shp \code{sf} polygons (or anything coercible via \code{st_as_sf}).
#' @param delta candidate-point spacing.
#' @param phi numeric vector of spatial scale parameters to profile; if
#'   \code{NULL}, a default grid from \code{sqrt(min area)} to \code{extent/10}.
#' @param method point method: 1 = SSI, 2 = uniform, 3 = regular grid.
#' @param weighted logical; population-weighted aggregation using \code{pop_shp}.
#' @param pop_shp population-density \code{SpatRaster} (needed if \code{weighted}).
#' @param kappa Matern smoothness for the spatial kernel (0.5 default).
#' @param par0 optional starting values \code{c(beta, sigma2, phi)}.
#' @param control.mcmc list from \code{\link{control_mcmc}}.
#' @param phi_method how the spatial scale is estimated: \code{"grid"} (profile
#'   over the supplied \code{phi} grid, the robust default) or \code{"direct"}
#'   (optimise \code{phi} continuously inside the MCML objective; exponential
#'   kernel only). See the package vignette/PDF on the double-integral derivation.
#' @param nugget logical; if \code{TRUE} (requires \code{phi_method = "direct"})
#'   add an unstructured region-level term, fitting covariance
#'   \eqn{\sigma^2(R(\phi)+\nu I)} and estimating the relative nugget
#'   \eqn{\nu=\tau^2/\sigma^2} with a standard error. Absorbs overdispersion not
#'   explained by the spatial structure.
#' @param reanchor number of re-anchoring passes: after fitting, the latent field
#'   is re-simulated at the current optimum and the model refit, which keeps the
#'   importance weights near-uniform (raises the MC effective sample size). 0
#'   (default) fits once; 2-3 is usually ample.
#' @param rho,giveup point-generation controls.
#' @param nthreads OpenMP threads for the correlation build.
#' @param messages logical; print optimiser progress.
#' @return an object of class \code{"SDALGCP2"}.
#' @seealso \code{\link{mcml_fit}}, \code{\link{precompute_corr}}, \code{\link{sda_points}}
#' @examples
#' \donttest{
#' library(sf)
#' ## ---- simulate a lattice of regions and aggregated counts ----
#' set.seed(1)
#' bound <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20)))
#' shp   <- st_sf(geometry = st_make_grid(bound, n = c(8, 8)))
#' N     <- nrow(shp)
#'
#' pts   <- sda_points(shp, delta = 1.2, method = 3)      # regular grid points
#' phi_grid <- seq(1, 5, length.out = 8)
#' corr  <- precompute_corr(pts, phi_grid)
#' Sig   <- 0.5 * corr$R[, , which.min(abs(phi_grid - 2.5))]
#' x1    <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
#' pop   <- round(runif(N, 500, 3000))
#' y     <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
#'                             as.numeric(t(chol(Sig)) %*% rnorm(N))))
#' dat   <- data.frame(y = y, x1 = x1, pop = pop)
#'
#' ## ---- fit ----
#' ctrl <- control_mcmc(n.sim = 6000, burnin = 1500, thin = 6, h = 1.65 / N^(1/6))
#' fit  <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.2,
#'                  phi = phi_grid, method = 3, control.mcmc = ctrl)
#' summary(fit)
#'
#' ## ---- predict ----
#' pred_d <- predict(fit, type = "discrete",   sampler = "mcmc",    control.mcmc = ctrl)
#' pred_c <- predict(fit, type = "continuous", sampler = "laplace", cellsize = 1,
#'                   control.mcmc = ctrl)
#' }
#' @export
SDALGCP2 <- function(formula, data, my_shp, delta, phi = NULL, method = 1L,
                     weighted = FALSE, pop_shp = NULL, kappa = 0.5,
                     par0 = NULL, control.mcmc = NULL, phi_method = c("grid", "direct"),
                     nugget = FALSE, reanchor = 0L, rho = 0.55, giveup = 1000L,
                     nthreads = 0L, messages = FALSE) {
  phi_method <- match.arg(phi_method)
  if (!inherits(formula, "formula")) stop("'formula' must be a formula.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!inherits(my_shp, "sf")) my_shp <- sf::st_as_sf(my_shp)

  if (is.null(phi)) {
    areas <- as.numeric(sf::st_area(my_shp))
    bb <- sf::st_bbox(my_shp)
    min_phi <- sqrt(min(areas))
    max_phi <- min(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) / 10
    phi <- seq(min_phi, max_phi, length.out = 20)
  }
  if (!is.null(par0)) phi <- sort(unique(c(phi, utils::tail(par0, 1))))

  pts  <- sda_points(my_shp, delta, method = method, weighted = weighted,
                     pop_shp = pop_shp, rho = rho, giveup = giveup)
  corr <- precompute_corr(pts, phi, kappa = kappa, nthreads = nthreads)
  fit  <- mcml_fit(formula, data, corr, par0 = par0, control.mcmc = control.mcmc,
                   phi_method = phi_method, nugget = nugget, reanchor = reanchor,
                   messages = messages)
  fit$call <- match.call()
  fit
}
