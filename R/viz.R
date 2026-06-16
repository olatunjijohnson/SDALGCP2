# Post-fit visualisation: the maps and summaries a disease-mapping analyst wants
# after fitting an SDA-LGCP model. ggplot2-based; works for both region-level
# (discrete) and continuous predictions.

.discrete_var <- function(type) switch(type,
  RR = "pMean_RR", SE_RR = "pSD_RR", ARR = "pMean_ARR", SE_ARR = "pSD_ARR",
  stop("variable must be one of 'RR','SE_RR','ARR','SE_ARR' for discrete maps."))

.var_label <- function(v) switch(v,
  RR = "Incidence\nrelative risk", SE_RR = "SD of\nincidence RR",
  ARR = "Covariate-adjusted\nrelative risk", SE_RR_c = "SD",
  SE_ARR = "SD of\nadjusted RR", RR_c = "Relative risk", v)

#' Map a fitted SDALGCP2 prediction
#'
#' Maps a relative-risk surface (and its uncertainty) from
#' \code{\link{predict.SDALGCP2}}. Region-level predictions are drawn as choropleths;
#' continuous predictions as a raster, optionally masked to a boundary.
#'
#' @param x an object of class \code{"SDALGCP2_pred"}.
#' @param variable for discrete: one of \code{"RR"}, \code{"SE_RR"}, \code{"ARR"},
#'   \code{"SE_ARR"}; for continuous: \code{"RR"} or \code{"SE_RR"}.
#' @param bound optional \code{sf} boundary; continuous surfaces are masked to it
#'   and its outline overlaid.
#' @param midpoint optional value to centre a diverging colour scale (e.g. 1 for
#'   relative risk); if \code{NULL} a sequential viridis scale is used.
#' @param title optional plot title.
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @method plot SDALGCP2_pred
#' @export
plot.SDALGCP2_pred <- function(x, variable = "RR", bound = NULL,
                               midpoint = NULL, title = NULL, ...) {
  lab <- .var_label(if (x$type == "continuous" && variable == "RR") "RR_c" else variable)

  if (x$type == "discrete") {
    shp <- x$my_shp
    if (is.null(shp)) stop("This prediction has no polygon geometry to map.")
    if (!inherits(shp, "sf")) shp <- sf::st_as_sf(shp)
    shp$.v <- shp[[.discrete_var(variable)]]
    p <- ggplot2::ggplot(shp) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$.v), color = "grey60", linewidth = 0.1)
  } else {
    df <- data.frame(x = x$pred.loc[, 1], y = x$pred.loc[, 2],
                     .v = if (variable == "SE_RR") x$RRsd else x$RRmean)
    if (!is.null(bound)) {
      if (!inherits(bound, "sf")) bound <- sf::st_as_sf(bound)
      pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(bound))
      df <- df[lengths(sf::st_intersects(pts, sf::st_union(bound))) > 0, ]
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$.v)) +
      ggplot2::geom_tile() + ggplot2::coord_equal()
    if (!is.null(bound))
      p <- p + ggplot2::geom_sf(data = bound, inherit.aes = FALSE,
                                fill = NA, color = "black")
  }

  scale <- if (!is.null(midpoint)) {
    ggplot2::scale_fill_gradient2(name = lab, midpoint = midpoint,
      low = "#2166AC", mid = "grey95", high = "#B2182B")
  } else {
    ggplot2::scale_fill_viridis_c(name = lab)
  }
  p + scale + ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = NULL, y = NULL)
}

#' Map exceedance probabilities P(relative risk > threshold)
#'
#' @param x an \code{"SDALGCP2_pred"} object.
#' @param threshold a single relative-risk threshold.
#' @param bound optional \code{sf} boundary (continuous only).
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @export
map_exceedance <- function(x, threshold = 1, bound = NULL, ...) {
  stopifnot(inherits(x, "SDALGCP2_pred"))
  ex <- as.numeric(exceedance(x, threshold))
  lab <- sprintf("P(RR > %g)", threshold)
  if (x$type == "discrete") {
    shp <- x$my_shp; if (!inherits(shp, "sf")) shp <- sf::st_as_sf(shp)
    shp$.e <- ex
    ggplot2::ggplot(shp) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$.e), color = "grey60", linewidth = 0.1) +
      ggplot2::scale_fill_viridis_c(name = lab, limits = c(0, 1), option = "magma") +
      ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = NULL)
  } else {
    df <- data.frame(x = x$pred.loc[, 1], y = x$pred.loc[, 2], .e = ex)
    if (!is.null(bound)) {
      if (!inherits(bound, "sf")) bound <- sf::st_as_sf(bound)
      pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(bound))
      df <- df[lengths(sf::st_intersects(pts, sf::st_union(bound))) > 0, ]
    }
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$.e)) +
      ggplot2::geom_tile() + ggplot2::coord_equal() +
      ggplot2::scale_fill_viridis_c(name = lab, limits = c(0, 1), option = "magma") +
      ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = NULL)
    if (!is.null(bound)) g <- g + ggplot2::geom_sf(data = bound, inherit.aes = FALSE,
                                                   fill = NA, color = "black")
    g
  }
}

#' Profile likelihood and confidence interval for the spatial scale phi
#'
#' Spline-smoothed profile deviance for \code{phi}, with the
#' \code{coverage}-level confidence interval where the deviance crosses the
#' chi-squared cutoff. Improves on the loess-based interval of \code{SDALGCP}.
#'
#' @param object a fitted \code{"SDALGCP2"} object.
#' @param coverage confidence level.
#' @param plot logical; draw the deviance curve.
#' @return invisibly, a list with the interval and the smoothed profile; a
#'   \code{ggplot} is drawn when \code{plot = TRUE}.
#' @export
phi_profile <- function(object, coverage = 0.95, plot = TRUE) {
  stopifnot(inherits(object, "SDALGCP2"))
  # Direct fits estimate phi continuously: report a Wald interval from the Hessian.
  if (identical(object$phi_method, "direct")) {
    se <- sqrt(object$cov["phi", "phi"])
    z <- stats::qnorm(1 - (1 - coverage) / 2)
    ci <- object$phi_opt + c(-1, 1) * z * se
    if (plot) message(sprintf("Direct phi estimate: %.3g, %d%% Wald CI [%.3g, %.3g]",
                              object$phi_opt, round(100 * coverage), ci[1], ci[2]))
    return(invisible(list(ci = ci, phi = object$phi_opt, se = se)))
  }
  phi <- object$all_para$phi; ll <- object$all_para$value
  if (length(unique(phi)) < 4) {
    grid <- phi; dev <- -2 * (ll - max(ll)); ci <- range(phi)
  } else {
    f <- stats::splinefun(phi, ll, method = "natural")
    grid <- seq(min(phi), max(phi), length.out = 1000)
    sm <- f(grid); dev <- -2 * (sm - max(sm))
    cut <- stats::qchisq(coverage, 1)
    # CI = the (contiguous) set of phi where the deviance is below the cutoff;
    # this handles boundary maxima where only one crossing exists.
    inside <- which(dev <= cut)
    ci <- if (length(inside)) grid[range(inside)] else range(phi)
  }
  if (plot) {
    df <- data.frame(phi = grid, deviance = dev)
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$phi, .data$deviance)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = stats::qchisq(coverage, 1),
                          linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = ci, linetype = "dotted", color = "blue") +
      ggplot2::geom_vline(xintercept = object$phi_opt, color = "darkgreen") +
      ggplot2::labs(title = "Profile deviance for spatial scale phi",
        subtitle = sprintf("phi-hat = %.3g,  %d%% CI [%.3g, %.3g]",
                           object$phi_opt, round(100 * coverage), ci[1], ci[2]),
        x = expression(phi), y = "Deviance") +
      ggplot2::theme_minimal()
    print(g)
  }
  invisible(list(ci = ci, phi = grid, deviance = dev))
}

#' Coefficient plot of fixed effects (and sigma^2) with confidence intervals
#'
#' @param object a fitted \code{"SDALGCP2"} object.
#' @param level confidence level.
#' @param intercept logical; include the intercept.
#' @return a \code{ggplot} object.
#' @export
coef_plot <- function(object, level = 0.95, intercept = FALSE) {
  ci <- confint(object, level = level)
  est <- .sda_estimates(object)
  nm <- rownames(object$cov)
  df <- data.frame(term = factor(nm, levels = rev(nm)),
                   est = est, lo = ci[, 1], hi = ci[, 2])
  if (!intercept) df <- df[df$term != "(Intercept)", ]
  ggplot2::ggplot(df, ggplot2::aes(.data$est, .data$term)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                            height = 0.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(title = "Parameter estimates", x = "Estimate", y = NULL) +
    ggplot2::theme_minimal()
}

#' @export
plot.SDALGCP2 <- function(x, ...) phi_profile(x, ...)
