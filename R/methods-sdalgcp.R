# glm-like methods for objects returned by sdalgcp().

.strip_sdalgcp <- function(x) { cl <- class(x); class(x) <- cl[cl != "sdalgcp"]; x }

#' @method print sdalgcp
#' @export
print.sdalgcp <- function(x, ...) {
  cat("Spatially discrete LGCP fit  [", x$mode, "]\n", sep = "")
  cat("Call:  "); print(x$call)
  N <- if (!is.null(x$N)) x$N else length(x$y) / max(1, x$T %||% 1)
  npts <- if (!is.null(x$S.coord)) sum(vapply(x$S.coord, function(z) nrow(as.matrix(z$xy)), 0L)) else NA_integer_
  cat(sprintf("\n%d regions", N))
  if (!is.null(x$T)) cat(sprintf(" x %d times", x$T))
  if (!is.na(npts)) cat(sprintf(",  %d candidate points  (delta = %.3g)", npts, x$delta))
  cat("\n\nCoefficients:\n")
  s <- summary(.strip_sdalgcp(x))
  stats::printCoefmat(s$coefficients, has.Pvalue = TRUE, P.values = TRUE, digits = 3, signif.stars = TRUE)
  cat(sprintf("\nSpatial scale  phi = %.4g", x$phi_opt))
  if (!is.null(x$nu_opt)) cat(sprintf("\nTemporal range  nu  = %.4g", x$nu_opt))
  if (!is.null(x$nugget) && isTRUE(x$nugget)) cat(sprintf("\nRelative nugget  nu = %.4g", x$nu_opt))
  cat(sprintf("\nLog-likelihood = %.2f\n", x$llike_val_opt))
  if (!is.null(s$mc)) cat(sprintf("Monte Carlo effective sample size: %.0f%% of %d draws\n",
                                  100 * s$mc$ESS_frac, s$mc$B))
  invisible(x)
}

#' @method summary sdalgcp
#' @export
summary.sdalgcp <- function(object, ...) summary(.strip_sdalgcp(object), ...)

#' Predict relative risk from an sdalgcp fit
#'
#' Returns a prediction object carrying, for every location, the posterior mean and
#' standard error of the relative risk \code{RR} (\eqn{\exp(\eta)=\exp(d'\beta+S)})
#' and the covariate-adjusted relative risk \code{ARR} (\eqn{\exp(S)}). Map it with
#' \code{plot()} and get hotspot probabilities with \code{\link{exceedance}}.
#'
#' @param object an \code{"sdalgcp"} fit.
#' @param type \code{"discrete"} (region level, default) or \code{"continuous"}
#'   (a grid surface). Ignored for spatio-temporal fits.
#' @param sampler \code{"mcmc"} (default) or \code{"laplace"}.
#' @param cellsize grid spacing for \code{type = "continuous"}.
#' @param ... passed to the underlying predictor.
#' @return an object of class \code{"SDALGCP2_pred"} (spatial) or
#'   \code{"SDALGCP2_ST_pred"} (spatio-temporal). For discrete spatial fits
#'   \code{$my_shp} is an \code{sf} with \code{RR_mean}, \code{RR_se},
#'   \code{ARR_mean}, \code{ARR_se} columns.
#' @method predict sdalgcp
#' @export
predict.sdalgcp <- function(object, type = c("discrete", "continuous"),
                            sampler = c("mcmc", "laplace"), cellsize = NULL, ...) {
  obj <- .strip_sdalgcp(object)
  if (!is.null(object$T)) return(stats::predict(obj, ...))   # spatio-temporal
  type <- match.arg(type); sampler <- match.arg(sampler)
  stats::predict(obj, type = type, sampler = sampler, cellsize = cellsize, ...)
}

#' Map an sdalgcp fit
#'
#' Predicts and maps a chosen quantity. Works for spatial fits (discrete or
#' continuous) and spatio-temporal fits (select a \code{time}).
#'
#' @param x an \code{"sdalgcp"} fit.
#' @param what one of \code{"RR"} (relative risk, default), \code{"ARR"}
#'   (covariate-adjusted relative risk), \code{"RR_se"}, \code{"ARR_se"} or
#'   \code{"exceedance"}.
#' @param type \code{"discrete"} (default) or \code{"continuous"} (spatial fits).
#' @param time for spatio-temporal fits, the time to map (default: first; use
#'   \code{NULL} to facet all times).
#' @param threshold threshold for \code{what = "exceedance"}.
#' @param which for exceedance: \code{"ARR"} (default) or \code{"RR"}.
#' @param cellsize grid spacing for \code{type = "continuous"}.
#' @param sampler \code{"mcmc"} (default) or \code{"laplace"}.
#' @param ... passed to the mapping layer.
#' @return a \code{ggplot} object.
#' @method plot sdalgcp
#' @export
plot.sdalgcp <- function(x, what = c("RR", "ARR", "RR_se", "ARR_se", "exceedance"),
                         type = c("discrete", "continuous"), time = NULL,
                         threshold = 1, which = c("ARR", "RR"), cellsize = NULL,
                         sampler = c("mcmc", "laplace"), ...) {
  what <- match.arg(what); which <- match.arg(which); sampler <- match.arg(sampler)
  obj <- .strip_sdalgcp(x)

  if (!is.null(x$T)) {                                        # spatio-temporal
    pr <- stats::predict(obj)
    if (is.null(time)) time <- pr$times[1]
    return(plot(pr, time = time, what = what, threshold = threshold, which = which))
  }

  type <- match.arg(type)
  pr <- stats::predict(obj, type = type, sampler = sampler, cellsize = cellsize)
  bound <- if (type == "continuous") x$data_sf else NULL
  if (what == "exceedance")
    return(map_exceedance(pr, threshold = threshold, which = which, bound = bound))
  plot(pr, variable = what, bound = bound)
}

#' @method confint sdalgcp
#' @export
confint.sdalgcp <- function(object, parm, level = 0.95, ...)
  confint(.strip_sdalgcp(object), parm = parm, level = level, ...)
