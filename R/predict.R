# Prediction for fitted SDA-LGCP models (bottleneck B5 uses cross_cov_cpp).
# Discrete: region-level incidence and covariate-adjusted relative risk.
# Continuous: change-of-support prediction on a grid.
# Two latent samplers: "mcmc" (MALA, accurate default) and "laplace" (Gaussian
# approximation N(mode, Sigma.tilde) — no MCMC, near-instant).

.matern_R <- function(d, phi, kappa = 0.5) {
  u <- d / phi
  if (kappa == 0.5) return(exp(-u))
  if (kappa == 1.5) { a <- sqrt(3) * u; return((1 + a) * exp(-a)) }
  if (kappa == 2.5) { a <- sqrt(5) * u; return((1 + a + a^2 / 3) * exp(-a)) }
  exp(-u)
}

.colSDs <- function(X) {
  n <- nrow(X)
  sqrt(pmax(0, (colMeans(X^2) - colMeans(X)^2) * n / (n - 1)))
}

# Draw latent linear-predictor samples (nsim x N) at the fitted optimum.
.draw_latent <- function(object, sampler = "mcmc", control.mcmc = NULL) {
  mu0 <- object$mu; Sigma0 <- object$Sigma_mat_opt
  y <- object$y; m <- object$m
  if (is.null(control.mcmc)) control.mcmc <- object$control.mcmc
  if (sampler == "laplace") {
    lap <- laplace_mode_poisson_cpp(y, m, mu0, Sigma0)
    n.keep <- (control.mcmc$n.sim - control.mcmc$burnin) %/% control.mcmc$thin
    L <- t(chol(lap$Sigma_tilde))
    Z <- matrix(stats::rnorm(n.keep * length(mu0)), n.keep, length(mu0))
    sweep(Z %*% t(L), 2, lap$mode, "+")
  } else {
    laplace_sampling(mu0, Sigma0, y, m, control.mcmc)$samples
  }
}

#' Predict relative risk from a fitted SDALGCP2 model
#'
#' @param object an object of class \code{"SDALGCP2"} from \code{\link{SDALGCP2}}
#'   or \code{\link{mcml_fit}}.
#' @param type \code{"discrete"} for region-level inference or \code{"continuous"}
#'   for a spatially continuous surface.
#' @param sampler \code{"mcmc"} (MALA, default) or \code{"laplace"} (fast Gaussian
#'   approximation, no MCMC).
#' @param cellsize grid spacing for continuous prediction (ignored if
#'   \code{pred.loc} supplied).
#' @param pred.loc optional data frame of prediction coordinates (\code{x},
#'   \code{y}) for continuous prediction.
#' @param control.mcmc optional MCMC controls; defaults to those used at fitting.
#' @param ... unused.
#' @return for \code{type = "discrete"}, an \code{sf} object augmented with
#'   posterior mean/SD of incidence relative risk (\code{pMean_RR}/\code{pSD_RR})
#'   and covariate-adjusted relative risk (\code{pMean_ARR}/\code{pSD_ARR}); for
#'   \code{type = "continuous"}, a list with the prediction grid and posterior
#'   summaries. Result carries class \code{"SDALGCP2_pred"}.
#' @method predict SDALGCP2
#' @export
predict.SDALGCP2 <- function(object, type = c("discrete", "continuous"),
                             sampler = c("mcmc", "laplace"), cellsize = NULL,
                             pred.loc = NULL, control.mcmc = NULL, ...) {
  type <- match.arg(type); sampler <- match.arg(sampler)
  mu0 <- object$mu
  S.sim <- .draw_latent(object, sampler, control.mcmc)

  if (type == "discrete") {
    shp <- attr(object, "my_shp")
    RR  <- exp(S.sim)
    ARR <- exp(sweep(S.sim, 2, mu0, "-"))
    if (!is.null(shp)) {
      shp$pMean_RR  <- colMeans(RR);  shp$pSD_RR  <- .colSDs(RR)
      shp$pMean_ARR <- colMeans(ARR); shp$pSD_ARR <- .colSDs(ARR)
    }
    out <- list(type = "discrete", my_shp = shp,
                pMean_RR = colMeans(RR), pSD_RR = .colSDs(RR),
                pMean_ARR = colMeans(ARR), pSD_ARR = .colSDs(ARR),
                S.draw = S.sim, mu = mu0)
    class(out) <- "SDALGCP2_pred"
    return(out)
  }

  ## continuous (change of support)
  shp <- attr(object, "my_shp")
  weighted <- isTRUE(attr(object, "weighted"))
  kappa <- if (!is.null(object$kappa)) object$kappa else 0.5
  sigma2 <- object$sigma2_opt; phi <- object$phi_opt
  S.coord <- object$S.coord
  if (is.null(pred.loc)) {
    if (is.null(cellsize)) stop("Provide 'cellsize' or 'pred.loc' for continuous prediction.")
    grid <- sf::st_make_grid(sf::st_union(shp), cellsize = cellsize,
                             what = "centers", square = TRUE)
    co <- sf::st_coordinates(grid)
    pred.loc <- data.frame(x = co[, 1], y = co[, 2])
  }
  pl <- as.matrix(pred.loc[, 1:2])

  coords <- lapply(S.coord, function(z) as.matrix(z$xy)[, 1:2, drop = FALSE])
  wts <- if (weighted) lapply(S.coord, function(z) as.numeric(z$weight)) else list()

  Sigma_xx <- sigma2 * .matern_R(as.matrix(stats::dist(pl)), phi, kappa)
  Sigma_xA <- sigma2 * cross_cov_cpp(pl, coords, wts, phi, kappa, weighted, 0L)
  invA <- solve(object$Sigma_mat_opt)
  Bmat <- Sigma_xA %*% invA
  predVar <- Sigma_xx - Bmat %*% t(Sigma_xA)
  predVar <- (predVar + t(predVar)) / 2
  K <- t(chol(predVar + diag(1e-8, nrow(predVar))))

  nsim <- nrow(S.sim); npred <- nrow(pl)
  Sx <- matrix(0, nsim, npred)
  for (i in seq_len(nsim)) {
    Sx[i, ] <- Bmat %*% (S.sim[i, ] - mu0) + K %*% stats::rnorm(npred)
  }
  out <- list(type = "continuous", pred.loc = pred.loc,
              RRmean = colMeans(exp(Sx)), RRsd = .colSDs(exp(Sx)),
              pred.draw = Sx, my_shp = shp)
  class(out) <- "SDALGCP2_pred"
  out
}

#' Exceedance probabilities P(relative risk > threshold)
#'
#' @param object an \code{"SDALGCP2_pred"} object from \code{\link{predict.SDALGCP2}}.
#' @param thresholds numeric vector of thresholds.
#' @return a matrix of exceedance probabilities (locations x thresholds).
#' @export
exceedance <- function(object, thresholds) {
  stopifnot(inherits(object, "SDALGCP2_pred"))
  draws <- if (object$type == "continuous") exp(object$pred.draw) else
    exp(sweep(object$S.draw, 2, object$mu, "-"))
  sapply(thresholds, function(t) apply(draws, 2, function(x) mean(x > t)))
}
