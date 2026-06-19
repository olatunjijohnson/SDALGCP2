# A single, intuitive entry point for fitting SDA-LGCP models, modelled on the
# familiarity of glm(): a formula, an sf data object, sensible defaults. It
# dispatches to the spatial, raster-covariate or spatio-temporal engine based on
# the inputs, so users do not need to know which low-level function to call.

#' Control settings for \code{\link{sdalgcp}}
#'
#' Bundles the technical knobs so that a default fit needs none of them.
#'
#' @param delta candidate-point spacing. If \code{NULL} (default) it is chosen
#'   automatically to place roughly \code{points_per_region} points in a typical
#'   region.
#' @param points_per_region target number of candidate points per region used to
#'   pick \code{delta} automatically.
#' @param point_method how candidate points are laid out: \code{"regular"}
#'   (deterministic grid, default), \code{"uniform"} or \code{"ssi"}.
#' @param scale how the spatial scale \eqn{\phi} is estimated: \code{"continuous"}
#'   (optimised directly, no grid -- the default) or \code{"grid"} (profiled over
#'   \code{phi}). Spatio-temporal fits always profile \eqn{\phi} on a grid.
#' @param phi optional \eqn{\phi} grid (only used when \code{scale = "grid"} or for
#'   spatio-temporal fits); chosen from the geometry if \code{NULL}.
#' @param kappa spatial Matern smoothness (\code{0.5}, \code{1.5} or \code{2.5}).
#' @param kappa_t temporal Matern smoothness (spatio-temporal fits).
#' @param nugget logical; add an unstructured region-level term (overdispersion).
#'   Requires \code{scale = "continuous"}.
#' @param confounding \code{"none"} (default) or \code{"restricted"}. With
#'   \code{"restricted"}, restricted spatial regression is used: the spatial random
#'   effect is constrained to the orthogonal complement of the fixed-effect design
#'   so it cannot absorb a spatially structured covariate (avoids spatial
#'   confounding / attenuation of \code{beta}). Spatial models only.
#' @param reanchor number of re-anchoring passes (re-simulate the latent field at
#'   the optimum and refit) for reliable variance estimates. Default \code{2}.
#' @param n_sim,burnin,thin MCMC length controls for the latent-field sampler.
#' @param tilt_spatial logical; for raster covariates, use the fully
#'   covariate-tilted correlation (see \code{\link{SDALGCP2_raster}}).
#' @param nthreads OpenMP threads for the correlation assembly (0 = default).
#' @return a list of control settings.
#' @seealso \code{\link{sdalgcp}}
#' @examples
#' ## defaults, then a faster grid-based fit with a nugget term
#' str(sdalgcp_control())
#' ctrl <- sdalgcp_control(scale = "grid", nugget = FALSE, n_sim = 4000,
#'                         burnin = 1000, thin = 6)
#' @export
sdalgcp_control <- function(delta = NULL, points_per_region = 16,
                            point_method = c("regular", "uniform", "ssi"),
                            scale = c("continuous", "grid"), phi = NULL,
                            kappa = 0.5, kappa_t = 0.5, nugget = FALSE,
                            confounding = c("none", "restricted"),
                            reanchor = 2L, n_sim = 10000L, burnin = 2000L,
                            thin = 8L, tilt_spatial = FALSE, nthreads = 0L) {
  point_method <- match.arg(point_method)
  scale <- match.arg(scale)
  confounding <- match.arg(confounding)
  stopifnot((n_sim - burnin) %% thin == 0)
  list(delta = delta, points_per_region = points_per_region,
       point_method = point_method, scale = scale, phi = phi, kappa = kappa,
       kappa_t = kappa_t, nugget = nugget, confounding = confounding,
       reanchor = as.integer(reanchor), n_sim = as.integer(n_sim),
       burnin = as.integer(burnin), thin = as.integer(thin),
       tilt_spatial = tilt_spatial, nthreads = nthreads)
}

# Auto candidate-point spacing: aim for ~points_per_region points in a typical
# (median-area) region. For points on a grid, n ~ area / delta^2.
.auto_delta <- function(shp, points_per_region = 16) {
  a <- stats::median(as.numeric(sf::st_area(shp)))
  as.numeric(sqrt(a / points_per_region))
}

#' Fit a spatially discrete LGCP model for aggregated counts
#'
#' The main user interface, designed to feel like \code{\link[stats]{glm}}: give a
#' formula and an \code{sf} data object and it does the rest. The same call covers
#' three settings, chosen from the arguments you supply:
#' \itemize{
#'   \item \strong{spatial} (default): \code{sdalgcp(y ~ x + offset(log(pop)), data)};
#'   \item \strong{raster covariates}: add \code{rasters =} a \code{SpatRaster}
#'     whose layers are named in the formula -- these enter on the intensity scale
#'     (see \code{\link{SDALGCP2_raster}});
#'   \item \strong{spatio-temporal}: add \code{time =} the name of a time column.
#' }
#'
#' @param formula a model formula, e.g. \code{cases ~ x1 + offset(log(pop))}.
#' @param data an \code{sf} object of polygons whose columns hold the response,
#'   covariates and offset (one row per region, or per region-time for
#'   spatio-temporal fits).
#' @param time optional name of a time column in \code{data}; if given, a
#'   spatio-temporal model is fitted (data must have one row per region and time).
#' @param rasters optional \code{terra::SpatRaster} of spatially continuous
#'   covariates (layers named in \code{formula}).
#' @param covariates optional named list of \code{sf} \strong{point} layers giving
#'   covariates observed on a different support (e.g. monitors); each is kriged to
#'   the candidate points and enters with a Berkson correction (see
#'   \code{\link{SDALGCP2_misaligned}}).
#' @param popden optional population-density \code{SpatRaster}; if supplied, the
#'   region aggregation is population-weighted.
#' @param control a \code{\link{sdalgcp_control}} list of settings (smart defaults).
#' @param verbose logical; print progress.
#' @return a fitted model object of class \code{c("sdalgcp", ...)} with
#'   \code{print}, \code{summary}, \code{confint}, \code{predict} and \code{plot}
#'   methods.
#' @seealso \code{\link{predict.sdalgcp}}, \code{\link{sdalgcp_control}},
#'   \code{\link{SDALGCP2}}, \code{\link{SDALGCP2_raster}}, \code{\link{SDALGCP2_ST}}
#' @examples
#' \donttest{
#' library(sf)
#' set.seed(1)
#' grid <- st_make_grid(st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))),
#'                      n = c(8, 8))
#' regions <- st_sf(geometry = grid)
#' regions$x1  <- as.numeric(scale(st_coordinates(st_centroid(regions))[, 1]))
#' regions$pop <- round(runif(nrow(regions), 500, 3000))
#' regions$cases <- rpois(nrow(regions), regions$pop * exp(-6 + 0.5 * regions$x1))
#'
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = regions)  # that's it
#' summary(fit)
#' rr <- predict(fit)          # an sf you can plot() directly
#' plot(fit)                   # default relative-risk map
#' }
#' @export
sdalgcp <- function(formula, data, time = NULL, rasters = NULL, covariates = NULL,
                    popden = NULL, control = sdalgcp_control(), verbose = FALSE) {
  if (!inherits(formula, "formula")) stop("'formula' must be a formula, e.g. cases ~ x + offset(log(pop)).")
  if (!inherits(data, "sf")) stop("'data' must be an 'sf' object (polygons with the model columns). Convert with sf::st_as_sf().")
  weighted <- !is.null(popden)
  pm <- c(ssi = 1L, uniform = 2L, regular = 3L)[[control$point_method]]
  ctrl_mcmc <- control_mcmc(n.sim = control$n_sim, burnin = control$burnin, thin = control$thin)

  ## ---- spatio-temporal ----
  if (!is.null(time)) {
    if (!time %in% names(data)) stop("'time' column '", time, "' not found in data.")
    df <- sf::st_drop_geometry(data)
    times <- sort(unique(df[[time]])); T <- length(times)
    if (nrow(data) %% T != 0) stop("Spatio-temporal data must have one row per region and time (nrow divisible by #times).")
    ord <- order(match(df[[time]], times))            # sort by time
    data <- data[ord, ]; df <- df[ord, ]
    shp <- data[df[[time]] == times[1], ]             # N region geometries (first time)
    delta <- control$delta %||% .auto_delta(shp, control$points_per_region)
    if (verbose) message(sprintf("Spatio-temporal fit: %d regions x %d times, delta = %.3g", nrow(shp), T, delta))
    fit <- SDALGCP2_ST(formula, df, shp, times = times, delta = delta, phi = control$phi,
                       kappa = control$kappa, kappa_t = control$kappa_t, method = pm,
                       weighted = weighted, pop_shp = popden, control.mcmc = ctrl_mcmc,
                       reanchor = control$reanchor, rasters = rasters, covariates = covariates,
                       confounding = control$confounding, messages = verbose)
    fit$mode <- paste0("spatio-temporal",
                       if (!is.null(rasters)) " (raster covariates)"
                       else if (!is.null(covariates)) " (misaligned covariates)"
                       else if (control$confounding == "restricted") " (restricted)" else "")
  } else if (!is.null(covariates)) {
    ## ---- covariates on a different support (kriged + Berkson) ----
    delta <- control$delta %||% .auto_delta(data, control$points_per_region)
    if (verbose) message(sprintf("Misaligned-covariate fit: %d regions, delta = %.3g", nrow(data), delta))
    fit <- SDALGCP2_misaligned(formula, data, delta = delta, covariates = covariates,
                               phi = control$phi, method = pm, weighted = weighted,
                               pop_shp = popden, control.mcmc = ctrl_mcmc, messages = verbose)
    fit$mode <- "spatial (misaligned covariates)"
  } else if (!is.null(rasters)) {
    ## ---- raster (intensity-scale) covariates ----
    df <- sf::st_drop_geometry(data)
    delta <- control$delta %||% .auto_delta(data, control$points_per_region)
    if (verbose) message(sprintf("Raster-covariate fit: %d regions, delta = %.3g", nrow(data), delta))
    fit <- SDALGCP2_raster(formula, df, data, delta = delta, rasters = rasters, phi = control$phi,
                           method = pm, weighted = weighted, pop_shp = popden,
                           kappa = control$kappa, tilt_spatial = control$tilt_spatial,
                           control.mcmc = ctrl_mcmc, messages = verbose)
    fit$mode <- "spatial (raster covariates)"
  } else {
    ## ---- spatial ----
    df <- sf::st_drop_geometry(data)
    delta <- control$delta %||% .auto_delta(data, control$points_per_region)
    if (verbose) message(sprintf("Spatial fit: %d regions, delta = %.3g, scale = %s", nrow(data), delta, control$scale))
    fit <- SDALGCP2(formula, df, data, delta = delta, phi = control$phi,
                    method = pm, weighted = weighted, pop_shp = popden, kappa = control$kappa,
                    control.mcmc = ctrl_mcmc, phi_method = if (control$scale == "continuous") "direct" else "grid",
                    nugget = control$nugget, confounding = control$confounding,
                    reanchor = control$reanchor, nthreads = control$nthreads, messages = verbose)
    fit$mode <- if (control$confounding == "restricted") "spatial (restricted)" else "spatial"
  }
  fit$delta <- delta; fit$formula <- formula; fit$call <- match.call()
  fit$data_sf <- if (is.null(time)) data else shp
  class(fit) <- c("sdalgcp", class(fit))
  fit
}

# null-coalescing helper (base R has %||% only from 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a
