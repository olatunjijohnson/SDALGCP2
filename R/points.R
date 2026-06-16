# Candidate-point generation for the SDA-LGCP discretisation (lean: sf + spatstat
# + terra only; drops the legacy splancs/sp paths of SDALGCP).
#
# For each region we place candidate points (regular grid, uniform random, or
# Simple Sequential Inhibition) and attach aggregation weights. Unweighted:
# w_ik = 1/n_i (the region average). Weighted: w_ik = rho(x_ik)/sum_k rho(x_ik),
# the Monte-Carlo / quadrature estimate of the population-weighted average
# T_i(S) = int_Ai S(x) rho(x) dx / int_Ai rho(x) dx, with rho read from a raster.

# Number of points from area and spacing (matches SDALGCP's packing heuristic).
.n_from_delta <- function(area, delta, rho = 0.55) {
  max(1L, round((rho * area * 4) / (pi * delta^2)))
}

# Generate the raw coordinates for a single polygon geometry (sfc of length 1).
.region_xy <- function(geom, delta, method, rho, giveup) {
  win <- spatstat.geom::as.owin(sf::st_make_valid(geom))
  if (method == 3L) {
    bb <- sf::st_bbox(geom)
    gx <- seq(bb["xmin"], bb["xmax"], by = delta)
    gy <- seq(bb["ymin"], bb["ymax"], by = delta)
    g  <- expand.grid(x = gx, y = gy)
    inside <- spatstat.geom::inside.owin(g$x, g$y, win)
    xy <- as.matrix(g[inside, , drop = FALSE])
    if (nrow(xy) == 0L) {              # fall back to centroid if grid too coarse
      ct <- sf::st_coordinates(sf::st_centroid(geom))[, 1:2, drop = FALSE]
      xy <- matrix(as.numeric(ct), 1, 2)
    }
  } else {
    area <- as.numeric(sf::st_area(geom))
    n <- .n_from_delta(area, delta, rho)
    pts <- if (method == 1L) {
      spatstat.random::rSSI(r = delta, n = n, win = win, giveup = giveup)
    } else {
      spatstat.random::runifpoint(n, win = win)
    }
    xy <- cbind(pts$x, pts$y)
  }
  colnames(xy) <- c("x", "y")
  xy
}

#' Generate candidate sampling points inside each region
#'
#' Lean replacement for \code{SDALGCP::SDALGCPpolygonpoints()}: for every polygon
#' feature in \code{my_shp} it produces candidate points and aggregation weights,
#' in the list format consumed by \code{\link{precompute_corr}}.
#'
#' @param my_shp an \code{sf} object of \code{POLYGON}/\code{MULTIPOLYGON} features.
#' @param delta point spacing (grid step / SSI inhibition distance).
#' @param method 1 = SSI (default), 2 = uniform random, 3 = regular grid.
#' @param weighted logical; if \code{TRUE}, weights are population density read
#'   from \code{pop_shp}, otherwise equal weights.
#' @param pop_shp a \code{terra::SpatRaster} of population density (required when
#'   \code{weighted = TRUE}).
#' @param rho packing density used to choose the number of points.
#' @param giveup SSI rejection limit.
#' @return a list of length \code{nrow(my_shp)}; each element has \code{xy} and
#'   \code{weight}. Carries \code{"weighted"} and \code{"my_shp"} attributes.
#' @export
sda_points <- function(my_shp, delta, method = 1L, weighted = FALSE,
                       pop_shp = NULL, rho = 0.55, giveup = 1000L) {
  if (!inherits(my_shp, "sf")) my_shp <- sf::st_as_sf(my_shp)
  if (weighted && is.null(pop_shp))
    stop("'pop_shp' (population raster) is required when weighted = TRUE.")
  if (weighted && !inherits(pop_shp, "SpatRaster")) pop_shp <- terra::rast(pop_shp)

  N <- nrow(my_shp)
  geoms <- sf::st_geometry(my_shp)
  crs <- sf::st_crs(my_shp)
  pb <- progress::progress_bar$new(
    format = "  generating points  :current/:total [:bar] :percent",
    total = N, width = 70, clear = FALSE)

  out <- vector("list", N)
  for (i in seq_len(N)) {
    xy <- .region_xy(geoms[i], delta, as.integer(method), rho, giveup)
    if (weighted) {
      v <- terra::extract(pop_shp,
             terra::vect(sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                         coords = c("x", "y"), crs = crs)), ID = FALSE)[[1]]
      v[is.na(v)] <- 0
      w <- if (sum(v) > 0) v / sum(v) else rep(1 / nrow(xy), nrow(xy))
    } else {
      w <- rep(1 / nrow(xy), nrow(xy))
    }
    out[[i]] <- list(xy = xy, weight = as.numeric(w))
    pb$tick()
  }
  attr(out, "weighted") <- weighted
  attr(out, "my_shp") <- my_shp
  out
}
