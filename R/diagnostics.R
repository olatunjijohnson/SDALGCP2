# Model-checking diagnostics for fitted SDA-LGCP models: observed-vs-fitted,
# Pearson residuals, and a residual Moran's I test for leftover spatial
# autocorrelation (computed without spdep, from queen contiguity).

# Row-standardised queen-contiguity Moran's I with a permutation p-value.
.moran <- function(z, shp, nsim = 999) {
  nb <- sf::st_touches(shp)               # neighbour index list
  n <- length(z); zc <- z - mean(z)
  num <- 0; S0 <- 0
  for (i in seq_len(n)) {
    j <- nb[[i]]
    if (length(j) == 0) next
    w <- 1 / length(j)                    # row-standardised
    num <- num + w * sum(zc[i] * zc[j])
    S0 <- S0 + w * length(j)
  }
  I <- (n / S0) * num / sum(zc^2)
  perm <- vapply(seq_len(nsim), function(b) {
    zp <- sample(zc); np <- 0
    for (i in seq_len(n)) {
      j <- nb[[i]]; if (length(j) == 0) next
      np <- np + (1 / length(j)) * sum(zp[i] * zp[j])
    }
    (n / S0) * np / sum(zp^2)
  }, numeric(1))
  list(I = I, expected = -1 / (n - 1),
       p_value = (1 + sum(perm >= I)) / (nsim + 1))
}

#' Posterior-predictive model checking for an SDALGCP2 fit
#'
#' Compares observed counts with fitted Poisson means, returns Pearson residuals,
#' and tests them for residual spatial autocorrelation with Moran's I. A
#' non-significant Moran's I indicates the spatial random effect has absorbed the
#' spatial structure.
#'
#' @param object a fitted \code{"SDALGCP2"} object.
#' @param pred a discrete prediction from \code{predict(object, "discrete")}; if
#'   \code{NULL} one is computed with the fitting MCMC controls.
#' @param nsim permutations for the Moran's I p-value.
#' @param plot logical; draw the observed-vs-fitted scatter.
#' @return invisibly, a list with \code{fitted}, \code{residuals} and \code{moran}.
#' @export
model_check <- function(object, pred = NULL, nsim = 999, plot = TRUE) {
  stopifnot(inherits(object, "SDALGCP2"))
  if (is.null(pred)) pred <- predict(object, type = "discrete")
  shp <- attr(object, "my_shp")
  if (is.null(shp)) stop("Model has no polygon geometry for spatial diagnostics.")
  if (!inherits(shp, "sf")) shp <- sf::st_as_sf(shp)

  y <- object$y
  fitted <- object$m * pred$relative_risk       # Poisson mean = pop * exp(eta)
  resid <- (y - fitted) / sqrt(fitted)
  mI <- .moran(resid, shp, nsim = nsim)

  if (plot) {
    df <- data.frame(fitted = fitted, observed = y)
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$fitted, .data$observed)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::labs(title = "Observed vs fitted counts",
        subtitle = sprintf("Residual Moran's I = %.3f (E = %.3f), p = %.3f",
                           mI$I, mI$expected, mI$p_value),
        x = "Fitted", y = "Observed") +
      ggplot2::theme_minimal()
    print(g)
  }
  invisible(list(fitted = fitted, residuals = resid, moran = mI))
}

#' Importance-sampling diagnostics for an MCML fit
#'
#' The MCML estimate reweights latent samples drawn at the anchor towards the
#' optimum. When the optimum is far from the anchor the weights become uneven and
#' the estimate unreliable. This reports the effective sample size of the
#' importance weights at the maximiser and a Monte Carlo standard error for the
#' maximised log-likelihood, \eqn{\mathrm{SE}\approx\sqrt{1/\mathrm{ESS}-1/B}}.
#'
#' @param object a fitted \code{"SDALGCP2"} object.
#' @param warn_frac warn if the ESS falls below this fraction of \eqn{B}.
#' @return invisibly, a list with \code{B}, \code{ESS}, \code{ESS_frac} and
#'   \code{se_loglik}.
#' @export
mc_diagnostics <- function(object, warn_frac = 0.1) {
  stopifnot(inherits(object, "SDALGCP2"))
  S <- object$S; D <- object$D; y <- object$y; m <- object$m
  n <- length(y); p <- ncol(D); B <- nrow(S)
  data_ll <- as.numeric(S %*% y) - as.numeric(exp(S) %*% m)

  # Use the anchor at which the *retained* latent field was actually drawn
  # (after re-anchoring this differs from the original par0).
  anchor <- if (!is.null(object$anchor)) object$anchor else object$par0
  beta0 <- anchor[1:p]; s20 <- anchor[p + 1]; phi0 <- anchor[p + 2]
  corr <- attr(object, "prematrix")
  if (identical(object$phi_method, "direct")) {
    pts <- attr(object, "S_coord"); wtd <- isTRUE(attr(object, "weighted"))
    cds <- lapply(pts, function(z) as.matrix(z$xy)[, 1:2, drop = FALSE])
    wts <- if (wtd) lapply(pts, function(z) as.numeric(z$weight)) else list()
    R0 <- corr_and_grad_cpp(cds, wts, phi0, object$kappa, wtd, 0L)$R
  } else {
    R0 <- corr$R[, , which.min(abs(corr$phi - phi0))]
  }
  c0 <- chol(R0); Den <- .mc_num_loglik(c(beta0, log(s20)), D, chol2inv(c0),
                                        2 * sum(log(diag(c0))), S, data_ll, n, p)$num

  Rb <- object$Sigma_mat_opt / object$sigma2_opt
  cb <- chol(Rb); num_opt <- .mc_num_loglik(c(object$beta_opt, log(object$sigma2_opt)),
                D, chol2inv(cb), 2 * sum(log(diag(cb))), S, data_ll, n, p)$num

  a <- num_opt - Den; a <- a - max(a)
  w <- exp(a); w <- w / sum(w)
  ess <- 1 / sum(w^2)
  out <- list(B = B, ESS = ess, ESS_frac = ess / B,
              se_loglik = sqrt(max(0, 1 / ess - 1 / B)))
  if (out$ESS_frac < warn_frac)
    warning(sprintf("Low importance-sampling ESS (%.0f of %d, %.1f%%): consider re-anchoring (iterate = TRUE) or a par0 closer to the optimum.",
                    ess, B, 100 * out$ESS_frac))
  invisible(out)
}

#' One-call panel of post-fit graphics
#'
#' Returns the maps and summaries an analyst usually wants after fitting:
#' relative-risk and uncertainty maps, an exceedance map, the coefficient plot
#' and the phi profile. The pieces are returned as a named list of \code{ggplot}
#' objects so they can be arranged or printed individually.
#'
#' @param object a fitted \code{"SDALGCP2"} object.
#' @param pred optional discrete prediction; computed if \code{NULL}.
#' @param threshold relative-risk threshold for the exceedance map.
#' @param ... passed to \code{\link{predict.SDALGCP2}} when \code{pred} is computed.
#' @return a named list of \code{ggplot} objects.
#' @export
report <- function(object, pred = NULL, threshold = 1.5, ...) {
  if (is.null(pred)) pred <- predict(object, type = "discrete", ...)
  list(
    relative_risk = plot(pred, variable = "adjusted_rr", title = "Covariate-adjusted relative risk"),
    uncertainty   = plot(pred, variable = "adjusted_rr_se", title = "Uncertainty (SE)"),
    exceedance    = map_exceedance(pred, threshold = threshold),
    coefficients  = coef_plot(object),
    phi_profile   = phi_profile(object, plot = FALSE)
  )
}
