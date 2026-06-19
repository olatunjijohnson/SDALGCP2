# Post-fit visualisation: the maps and summaries a disease-mapping analyst wants
# after fitting an SDA-LGCP model. ggplot2-based; works for both region-level
# (discrete) and continuous predictions.

.pred_columns <- c("relative_risk", "adjusted_rr", "relative_risk_se", "adjusted_rr_se")

.pred_values <- function(x, variable) {
  v <- x[[variable]]
  if (is.null(v))
    stop("variable must be one of ", paste0("'", .pred_columns, "'", collapse = ", "), ".")
  as.numeric(v)
}

.var_label <- function(v) switch(v,
  relative_risk    = "Relative risk", adjusted_rr      = "Covariate-adjusted\nrelative risk",
  relative_risk_se = "SE\n(relative risk)", adjusted_rr_se = "SE\n(adjusted RR)", v)

# Short, human-readable label for the exceedance quantity (column tokens are ugly).
.which_label <- function(w) switch(w,
  adjusted_rr = "adjusted RR", relative_risk = "RR", w)

#' Map a fitted SDALGCP2 prediction
#'
#' Maps any of the four predicted quantities from \code{\link{predict.SDALGCP2}}
#' -- the relative risk \code{"relative_risk"}, the covariate-adjusted relative
#' risk \code{"adjusted_rr"}, or their standard errors
#' \code{"relative_risk_se"}/\code{"adjusted_rr_se"} -- for either discrete
#' (choropleth) or continuous (raster) predictions.
#'
#' @param x an object of class \code{"SDALGCP2_pred"}.
#' @param variable one of \code{"relative_risk"}, \code{"adjusted_rr"},
#'   \code{"relative_risk_se"}, \code{"adjusted_rr_se"}.
#' @param bound optional \code{sf} boundary; continuous surfaces are masked to it
#'   and its outline overlaid.
#' @param midpoint optional value to centre a diverging colour scale (defaults to 1
#'   for the relative-risk columns, none for the standard errors).
#' @param title optional plot title.
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
#'                control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
#'                                          reanchor = 0))
#' pr <- predict(fit, type = "discrete")
#' plot(pr, variable = "relative_risk")          # choropleth of relative risk
#' plot(pr, variable = "adjusted_rr_se")         # its uncertainty
#' }
#' @method plot SDALGCP2_pred
#' @export
plot.SDALGCP2_pred <- function(x, variable = c("relative_risk", "adjusted_rr",
                                               "relative_risk_se", "adjusted_rr_se"),
                               bound = NULL, midpoint = NULL, title = NULL, ...) {
  variable <- match.arg(variable)
  lab <- .var_label(variable); vals <- .pred_values(x, variable)
  if (is.null(midpoint) && variable %in% c("relative_risk", "adjusted_rr")) midpoint <- 1

  if (identical(attr(x, "pred_type"), "discrete")) {
    if (!inherits(x, "sf")) stop("This prediction has no polygon geometry to map.")
    shp <- x
    shp$fillvalue <- vals
    p <- ggplot2::ggplot(shp) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$fillvalue), color = "grey60", linewidth = 0.1)
  } else {
    pl <- attr(x, "pred_loc")
    df <- data.frame(x = pl[, 1], y = pl[, 2], fillvalue = vals)
    if (!is.null(bound)) {
      if (!inherits(bound, "sf")) bound <- sf::st_as_sf(bound)
      pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(bound))
      df <- df[lengths(sf::st_intersects(pts, sf::st_union(bound))) > 0, ]
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$fillvalue)) +
      ggplot2::geom_tile() + ggplot2::coord_equal()
    if (!is.null(bound))
      p <- p + ggplot2::geom_sf(data = sf::st_union(bound), inherit.aes = FALSE,
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

#' Map exceedance probabilities P(risk > threshold)
#'
#' @param x an \code{"SDALGCP2_pred"} object.
#' @param threshold a single relative-risk threshold.
#' @param which \code{"adjusted_rr"} (covariate-adjusted, default) or
#'   \code{"relative_risk"}.
#' @param bound optional \code{sf} boundary (continuous only).
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @seealso \code{\link{exceedance}} for the underlying probabilities.
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
#'                control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
#'                                          reanchor = 0))
#' pr <- predict(fit, type = "discrete")
#' map_exceedance(pr, threshold = 1.5)           # P(adjusted RR > 1.5)
#' }
#' @export
map_exceedance <- function(x, threshold = 1, which = c("adjusted_rr", "relative_risk"),
                           bound = NULL, ...) {
  stopifnot(inherits(x, "SDALGCP2_pred"))
  which <- match.arg(which)
  ex <- as.numeric(exceedance(x, threshold, which = which))
  lab <- sprintf("P(%s > %g)", .which_label(which), threshold)
  if (identical(attr(x, "pred_type"), "discrete")) {
    if (!inherits(x, "sf")) stop("This prediction has no polygon geometry to map.")
    shp <- x; shp$fillvalue <- ex
    ggplot2::ggplot(shp) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$fillvalue), color = "grey60", linewidth = 0.1) +
      ggplot2::scale_fill_viridis_c(name = lab, limits = c(0, 1), option = "magma") +
      ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = NULL)
  } else {
    pl <- attr(x, "pred_loc")
    df <- data.frame(x = pl[, 1], y = pl[, 2], fillvalue = ex)
    if (!is.null(bound)) {
      if (!inherits(bound, "sf")) bound <- sf::st_as_sf(bound)
      pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(bound))
      df <- df[lengths(sf::st_intersects(pts, sf::st_union(bound))) > 0, ]
    }
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$fillvalue)) +
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
#' chi-squared cutoff.
#'
#' @param object a fitted \code{"SDALGCP2"} object.
#' @param coverage confidence level.
#' @param plot logical; draw the deviance curve.
#' @return invisibly, a list with the interval and the smoothed profile; a
#'   \code{ggplot} is drawn when \code{plot = TRUE}.
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' ## profile phi on a grid (scale = "grid") so there is a deviance curve to draw
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
#'                control = sdalgcp_control(scale = "grid", n_sim = 2000,
#'                                          burnin = 500, thin = 5, reanchor = 0))
#' phi_profile(fit)
#' }
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
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = sdalgcp_data,
#'                control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
#'                                          reanchor = 0))
#' coef_plot(fit)
#' }
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

#' Plot an SDALGCP2 fit (the phi profile deviance)
#' @param x an \code{"SDALGCP2"} object.
#' @param ... passed to \code{\link{phi_profile}}.
#' @return invisibly, the profile (see \code{\link{phi_profile}}).
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' fit <- SDALGCP2(cases ~ x1 + offset(log(pop)),
#'                 sf::st_drop_geometry(sdalgcp_data), sdalgcp_data, delta = 1.2,
#'                 control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
#' plot(fit)   # profile deviance for the spatial scale phi
#' }
#' @method plot SDALGCP2
#' @export
plot.SDALGCP2 <- function(x, ...) phi_profile(x, ...)
