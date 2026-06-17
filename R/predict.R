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
#' @return an object of class \code{"SDALGCP2_pred"}. For both \code{type}s it
#'   carries, for every location (region or grid cell), the posterior mean and
#'   standard error of two quantities:
#'   \describe{
#'     \item{\code{RR} (relative risk)}{\eqn{\exp(\eta)=\exp(d'\beta+S)} -- the full
#'       relative risk, including the covariate effect.}
#'     \item{\code{ARR} (covariate-adjusted relative risk)}{\eqn{\exp(S)} -- the
#'       residual spatial relative risk after adjusting for covariates.}
#'   }
#'   stored as \code{RR_mean}/\code{RR_se}/\code{ARR_mean}/\code{ARR_se} (and, for
#'   discrete fits, as columns of the returned \code{sf}). Posterior draws are kept
#'   so that \code{\link{exceedance}} can be computed for either quantity.
#' @method predict SDALGCP2
#' @export
predict.SDALGCP2 <- function(object, type = c("discrete", "continuous"),
                             sampler = c("mcmc", "laplace"), cellsize = NULL,
                             pred.loc = NULL, control.mcmc = NULL, ...) {
  type <- match.arg(type); sampler <- match.arg(sampler)
  mu0 <- object$mu
  S.sim <- .draw_latent(object, sampler, control.mcmc)     # eta draws (nsim x N)
  shp <- attr(object, "my_shp")
  if (!is.null(shp) && !inherits(shp, "sf")) shp <- sf::st_as_sf(shp)

  if (type == "discrete") {
    RR  <- exp(S.sim)                                      # exp(eta) = relative risk
    ARR <- exp(sweep(S.sim, 2, mu0, "-"))                  # exp(S)   = covariate-adjusted
    out <- list(type = "discrete", my_shp = shp, mu = mu0, eta.draw = S.sim,
                RR_mean = colMeans(RR), RR_se = .colSDs(RR),
                ARR_mean = colMeans(ARR), ARR_se = .colSDs(ARR))
    if (!is.null(shp)) {
      shp$RR_mean <- out$RR_mean; shp$RR_se <- out$RR_se
      shp$ARR_mean <- out$ARR_mean; shp$ARR_se <- out$ARR_se
      out$my_shp <- shp
    }
    class(out) <- "SDALGCP2_pred"
    return(out)
  }

  ## continuous (change of support): predict the spatial field S(x) on a grid
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

  # region linear predictor mu(x) at each grid point (for the incidence RR)
  mu_x <- rep(NA_real_, npred)
  if (!is.null(shp)) {
    psf <- sf::st_as_sf(data.frame(x = pl[, 1], y = pl[, 2]), coords = c("x", "y"),
                        crs = sf::st_crs(shp))
    reg <- as.integer(suppressMessages(sf::st_intersects(psf, shp)))
    inside <- !is.na(reg); mu_x[inside] <- mu0[reg[inside]]
    if (any(!inside)) mu_x[!inside] <- mu0[sf::st_nearest_feature(psf[!inside, ], shp)]
  } else mu_x <- rep(mean(mu0), npred)

  ARR <- exp(Sx)                                          # exp(S(x))
  RR  <- exp(sweep(Sx, 2, mu_x, "+"))                     # exp(mu(x) + S(x))
  out <- list(type = "continuous", pred.loc = pred.loc, my_shp = shp,
              mu_x = mu_x, pred.draw = Sx,
              ARR_mean = colMeans(ARR), ARR_se = .colSDs(ARR),
              RR_mean = colMeans(RR), RR_se = .colSDs(RR))
  class(out) <- "SDALGCP2_pred"
  out
}

#' Exceedance probabilities P(risk > threshold)
#'
#' @param object an \code{"SDALGCP2_pred"} object from \code{\link{predict.SDALGCP2}}.
#' @param thresholds numeric vector of thresholds.
#' @param which which quantity: \code{"ARR"} (covariate-adjusted relative risk,
#'   default) or \code{"RR"} (relative risk).
#' @return a matrix of exceedance probabilities (locations x thresholds).
#' @export
exceedance <- function(object, thresholds, which = c("ARR", "RR")) {
  stopifnot(inherits(object, "SDALGCP2_pred"))
  which <- match.arg(which)
  draws <- if (object$type == "continuous") {
    if (which == "ARR") exp(object$pred.draw) else exp(sweep(object$pred.draw, 2, object$mu_x, "+"))
  } else {
    if (which == "ARR") exp(sweep(object$eta.draw, 2, object$mu, "-")) else exp(object$eta.draw)
  }
  sapply(thresholds, function(t) apply(draws, 2, function(x) mean(x > t)))
}
