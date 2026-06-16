#' Print an SDALGCP2 fit
#' @param x an \code{"SDALGCP2"} object.
#' @param ... unused.
#' @return \code{x}, invisibly.
#' @method print SDALGCP2
#' @export
print.SDALGCP2 <- function(x, ...) {
  cat("SDA-LGCP fit (SDALGCP2)\n")
  cat("Call: "); print(x$call)
  cf <- if (!is.null(x$estimates)) x$estimates else stats::setNames(
    c(x$beta_opt, x$sigma2_opt), rownames(x$cov))
  cat("\nCoefficients:\n"); print(round(cf, 4))
  cat(sprintf("Spatial scale phi: %g (%s)\n", x$phi_opt,
              if (!is.null(x$phi_method)) x$phi_method else "grid"))
  cat(sprintf("Log-likelihood:   %g\n", x$llike_val_opt))
  invisible(x)
}

# Coefficient vector aligned to the rows/cols of the covariance matrix.
.sda_estimates <- function(object) {
  if (!is.null(object$estimates)) return(object$estimates)
  stats::setNames(c(object$beta_opt, object$sigma2_opt), rownames(object$cov))
}

#' Summary of an SDALGCP2 fit
#'
#' @param object an object of class \code{"SDALGCP2"} from \code{\link{mcml_fit}}.
#' @param ... unused.
#' @return an object of class \code{"summary.SDALGCP2"} with a coefficient table.
#' @method summary SDALGCP2
#' @export
summary.SDALGCP2 <- function(object, ...) {
  est <- .sda_estimates(object)
  se  <- sqrt(diag(object$cov))
  z   <- est / se
  tab <- cbind(Estimate = est, Std.Err = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(tab) <- rownames(object$cov)
  mc <- tryCatch(suppressWarnings(mc_diagnostics(object)), error = function(e) NULL)
  out <- list(coefficients = tab, phi = object$phi_opt,
              loglik = object$llike_val_opt, mc = mc, call = object$call)
  class(out) <- "summary.SDALGCP2"
  out
}

#' Print a summary of an SDALGCP2 fit
#' @param x a \code{"summary.SDALGCP2"} object.
#' @param ... unused.
#' @return \code{x}, invisibly.
#' @method print summary.SDALGCP2
#' @export
print.summary.SDALGCP2 <- function(x, ...) {
  cat("Call: "); print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, has.Pvalue = TRUE, P.values = TRUE, digits = 3)
  cat(sprintf("\nSpatial scale phi: %g\nLog-likelihood: %g\n", x$phi, x$loglik))
  if (!is.null(x$mc))
    cat(sprintf("MC importance-sampling ESS: %.0f / %d (%.0f%%);  log-lik MC SE: %.3g\n",
                x$mc$ESS, x$mc$B, 100 * x$mc$ESS_frac, x$mc$se_loglik))
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
  est <- .sda_estimates(object)
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
