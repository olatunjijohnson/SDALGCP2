#' @useDynLib SDALGCP2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Precompute aggregated region-level correlation matrices
#'
#' Builds the \eqn{N \times N \times} \code{length(phi)} array of region-level
#' correlations used by the SDA-LGCP model, where
#' \deqn{R(\phi)_{ij} = \sum_{k,l} w_{ik} w_{jl}\, C(\lVert x_{ik}-x_{jl}\rVert; \phi, \kappa)}
#' (population-weighted) or the unweighted mean over candidate-point pairs. The
#' heavy reduction runs in C++ (OpenMP-parallel over region pairs); this is the
#' modern, fast replacement for \code{SDALGCP::precomputeCorrMatrix()}.
#'
#' @param points a list of length \eqn{N}; each element holds \code{$xy} (an
#'   \eqn{n_i \times 2} matrix of candidate-point coordinates) and, when
#'   weighted, \code{$weight} (a length-\eqn{n_i} vector summing to 1). The
#'   \code{"weighted"} and \code{"my_shp"} attributes produced by the
#'   point-generation step are honoured and carried through.
#' @param phi numeric vector of spatial scale parameters.
#' @param kappa Matern smoothness; \code{0.5} (exponential, default), \code{1.5}
#'   or \code{2.5} use closed forms in C++.
#' @param weighted logical; if \code{NULL} (default) it is taken from
#'   \code{attr(points, "weighted")}.
#' @param nthreads number of OpenMP threads; \code{0} (default) uses the OpenMP
#'   runtime default.
#'
#' @return a list with \code{R} (the correlation array) and \code{phi}, carrying
#'   \code{weighted}, \code{my_shp} and \code{S_coord} attributes on \code{R}.
#' @seealso \code{\link{sda_points}}, \code{\link{mcml_fit}}
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' pts <- sda_points(sdalgcp_data, delta = 1.2, method = 3)
#' cc  <- precompute_corr(pts, phi = c(2, 4, 6))
#' dim(cc$R)            # N x N x length(phi)
#' }
#' @export
precompute_corr <- function(points, phi, kappa = 0.5, weighted = NULL,
                            nthreads = 0L) {
  if (is.null(weighted)) weighted <- isTRUE(attr(points, "weighted"))
  phi <- as.numeric(phi)
  if (any(phi <= 0)) stop("'phi' values must be positive.")

  coords <- lapply(points, function(z) {
    xy <- as.matrix(z$xy)
    storage.mode(xy) <- "double"
    xy[, 1:2, drop = FALSE]
  })
  wts <- if (weighted) {
    lapply(points, function(z) as.numeric(z$weight))
  } else {
    list()
  }

  R <- corr_aggregate_cpp(coords, wts, phi, kappa, weighted, as.integer(nthreads))
  attr(R, "weighted") <- weighted
  attr(R, "my_shp")  <- attr(points, "my_shp")
  attr(R, "S_coord") <- points
  attr(R, "kappa")   <- kappa
  list(R = R, phi = phi, kappa = kappa)
}

#' Reference (pure-R) aggregated correlation builder
#'
#' Slow but dependency-light implementation kept for correctness testing and
#' benchmarking against \code{\link{precompute_corr}}. Computes the exponential
#' (\code{kappa = 0.5}) aggregated correlation only.
#'
#' @inheritParams precompute_corr
#' @return the \eqn{N \times N \times} \code{length(phi)} correlation array.
#' @keywords internal
precompute_corr_ref <- function(points, phi, weighted = NULL) {
  if (is.null(weighted)) weighted <- isTRUE(attr(points, "weighted"))
  N <- length(points)
  P <- length(phi)
  R <- array(0, dim = c(N, N, P))
  for (i in seq_len(N)) {
    xi <- as.matrix(points[[i]]$xy)[, 1:2, drop = FALSE]
    wi <- if (weighted) points[[i]]$weight else NULL
    for (j in i:N) {
      xj <- as.matrix(points[[j]]$xy)[, 1:2, drop = FALSE]
      wj <- if (weighted) points[[j]]$weight else NULL
      D <- sqrt(outer(xi[, 1], xj[, 1], "-")^2 + outer(xi[, 2], xj[, 2], "-")^2)
      for (k in seq_len(P)) {
        K <- exp(-D / phi[k])
        val <- if (weighted) sum(outer(wi, wj) * K) else mean(K)
        R[i, j, k] <- R[j, i, k] <- val
      }
    }
  }
  R
}
