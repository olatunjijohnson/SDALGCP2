#' @export
print.SDALGCP2 <- function(x, ...) {
  cat("SDA-LGCP fit (SDALGCP2)\n")
  cat("Call: "); print(x$call)
  cf <- c(x$beta_opt, x$sigma2_opt)
  names(cf) <- rownames(x$cov)
  cat("\nCoefficients:\n"); print(round(cf, 4))
  cat(sprintf("Spatial scale phi: %g\n", x$phi_opt))
  cat(sprintf("Log-likelihood:   %g\n", x$llike_val_opt))
  invisible(x)
}

#' Summary of an SDALGCP2 fit
#'
#' @param object an object of class \code{"SDALGCP2"} from \code{\link{mcml_fit}}.
#' @param ... unused.
#' @return an object of class \code{"summary.SDALGCP2"} with a coefficient table.
#' @export
summary.SDALGCP2 <- function(object, ...) {
  est <- c(object$beta_opt, object$sigma2_opt)
  se  <- sqrt(diag(object$cov))
  z   <- est / se
  tab <- cbind(Estimate = est, Std.Err = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(tab) <- rownames(object$cov)
  out <- list(coefficients = tab, phi = object$phi_opt,
              loglik = object$llike_val_opt, call = object$call)
  class(out) <- "summary.SDALGCP2"
  out
}

#' @export
print.summary.SDALGCP2 <- function(x, ...) {
  cat("Call: "); print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, has.Pvalue = TRUE, P.values = TRUE, digits = 3)
  cat(sprintf("\nSpatial scale phi: %g\nLog-likelihood: %g\n", x$phi, x$loglik))
  cat("Note: sigma^2 is the variance of the latent Gaussian process.\n")
  invisible(x)
}

#' Wald confidence intervals for an SDALGCP2 fit
#'
#' @param object an object of class \code{"SDALGCP2"}.
#' @param parm parameters to report (names or indices); default all.
#' @param level confidence level.
#' @param ... unused.
#' @return a matrix of lower/upper confidence limits.
#' @export
confint.SDALGCP2 <- function(object, parm, level = 0.95, ...) {
  est <- c(object$beta_opt, object$sigma2_opt)
  nm  <- rownames(object$cov)
  names(est) <- nm
  se <- sqrt(diag(object$cov))
  if (missing(parm)) parm <- nm else if (is.numeric(parm)) parm <- nm[parm]
  a <- (1 - level) / 2
  q <- stats::qnorm(c(a, 1 - a))
  ci <- est[parm] + outer(se[parm], q)
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
  ci
}
