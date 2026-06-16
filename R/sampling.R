#' MCMC control settings for the MALA sampler
#'
#' @param n.sim total number of iterations.
#' @param burnin burn-in iterations to discard.
#' @param thin thinning interval; \code{(n.sim - burnin)} must be a multiple.
#' @param h initial Langevin step size; if missing, \code{1.65 / d^(1/6)} is used.
#' @param c1.h,c2.h step-size adaptation constants.
#' @return a named list consumed by \code{\link{laplace_sampling}} / the fit.
#' @export
control_mcmc <- function(n.sim = 10000, burnin = 2000, thin = 8,
                         h = NULL, c1.h = 0.01, c2.h = 1e-4) {
  if (is.null(h)) h <- Inf
  stopifnot(n.sim > burnin, thin > 0, (n.sim - burnin) %% thin == 0,
            h >= 0, c1.h >= 0, c2.h >= 0, c2.h <= 1)
  list(n.sim = n.sim, burnin = burnin, thin = thin, h = h, c1.h = c1.h, c2.h = c2.h)
}

#' Sample the latent field [S | Y] (Poisson, non-nested) via C++ MALA
#'
#' Fast drop-in for \code{SDALGCP::Laplace.sampling()} for the Poisson,
#' non-nested case: the Laplace mode (Newton) and the adaptive MALA loop both run
#' in C++. Given the same mode/covariance and seed it reproduces the original R
#' sampler bit-for-bit.
#'
#' @param mu prior mean vector.
#' @param Sigma prior covariance matrix.
#' @param y count vector.
#' @param units.m offset vector.
#' @param control.mcmc list from \code{\link{control_mcmc}}.
#' @return list with \code{samples} (kept x n matrix) and \code{h} (step sizes).
#' @export
laplace_sampling <- function(mu, Sigma, y, units.m, control.mcmc) {
  mu <- as.numeric(mu); y <- as.numeric(y); units.m <- as.numeric(units.m)
  lap <- laplace_mode_poisson_cpp(y, units.m, mu, Sigma)
  out <- mala_poisson_cpp(y, units.m, mu, Sigma, lap$mode, lap$Sigma_tilde,
                          control.mcmc$n.sim, control.mcmc$burnin, control.mcmc$thin,
                          control.mcmc$h, control.mcmc$c1.h, control.mcmc$c2.h)
  out
}
