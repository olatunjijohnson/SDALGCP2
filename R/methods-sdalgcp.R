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
#' Returns the fitted region-level relative risk as an \code{sf} object (for
#' spatial fits) so it can be mapped directly, or a long data frame for
#' spatio-temporal fits.
#'
#' @param object an \code{"sdalgcp"} fit.
#' @param type \code{"risk"} for covariate-adjusted relative risk
#'   \eqn{\exp(S)} (default), \code{"incidence"} for \eqn{\exp(\mu+S)}, or
#'   \code{"exceedance"} for \eqn{P(\mathrm{risk} > \mathrm{threshold})}.
#' @param threshold threshold for \code{type = "exceedance"}.
#' @param ... passed to the underlying predictor.
#' @return for spatial fits, the model's \code{sf} augmented with
#'   \code{relative_risk}, \code{relative_risk_se} (and \code{incidence},
#'   \code{exceedance} as requested); for spatio-temporal fits, a list with
#'   region-by-time matrices and a long table.
#' @method predict sdalgcp
#' @export
predict.sdalgcp <- function(object, type = c("risk", "incidence", "exceedance"),
                            threshold = 1, ...) {
  type <- match.arg(type)
  obj <- .strip_sdalgcp(object)
  if (!is.null(object$T)) return(stats::predict(obj, ...))  # spatio-temporal

  pr <- stats::predict(obj, type = "discrete", ...)
  shp <- object$data_sf
  if (!inherits(shp, "sf")) shp <- sf::st_as_sf(shp)
  shp$relative_risk    <- pr$pMean_ARR
  shp$relative_risk_se <- pr$pSD_ARR
  shp$incidence        <- pr$pMean_RR
  if (type == "exceedance") shp$exceedance <- as.numeric(exceedance(pr, threshold))
  shp
}

#' Map an sdalgcp fit
#'
#' Default visualisation: a choropleth of the covariate-adjusted relative risk
#' (spatial fits). Equivalent to \code{plot(predict(object), ...)}.
#'
#' @param x an \code{"sdalgcp"} fit.
#' @param type \code{"risk"} (default), \code{"incidence"}, \code{"risk_se"} or
#'   \code{"exceedance"}.
#' @param threshold threshold for \code{type = "exceedance"}.
#' @param ... passed to the mapping layer.
#' @return a \code{ggplot} object.
#' @method plot sdalgcp
#' @export
plot.sdalgcp <- function(x, type = c("risk", "incidence", "risk_se", "exceedance"),
                         threshold = 1, ...) {
  type <- match.arg(type)
  if (!is.null(x$T)) stop("Plotting spatio-temporal fits directly is not yet supported; use predict() and map per time slice.")
  obj <- .strip_sdalgcp(x)
  pr <- stats::predict(obj, type = "discrete")   # pr$my_shp already carries geometry + risk columns
  if (type == "exceedance") return(map_exceedance(pr, threshold = threshold, ...))
  var <- switch(type, risk = "ARR", incidence = "RR", risk_se = "SE_ARR")
  plot(pr, variable = var, midpoint = if (type %in% c("risk", "incidence")) 1 else NULL,
       title = switch(type, risk = "Covariate-adjusted relative risk",
                      incidence = "Incidence relative risk", risk_se = "Relative-risk SD"), ...)
}

#' @method confint sdalgcp
#' @export
confint.sdalgcp <- function(object, parm, level = 0.95, ...)
  confint(.strip_sdalgcp(object), parm = parm, level = level, ...)
