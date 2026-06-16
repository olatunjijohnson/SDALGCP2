# Spatio-temporal SDA-LGCP with a separable, KRONECKER-FREE likelihood (task #13b).
#
# Latent field over N regions x T times has separable covariance
#   Cov(vec M) = sigma2 (R_t(nu) %x% R_s(phi)),   M is N x T (regions x times).
# The original SDALGCP materialised the (NT)x(NT) Kronecker product and its
# inverse/log-det on every likelihood evaluation. We never form it, using
#   x' (R_t %x% R_s)^{-1} x = tr( R_s^{-1} M R_t^{-1} M' ),
#   log|R_t %x% R_s|        = N log|R_t| + T log|R_s|
# (both verified against the brute-force Kronecker computation), turning (NT)^3
# work into N^3 + T^3. The temporal correlation is Matern with range nu and
# smoothness kappa_t (estimated continuously); the spatial scale phi is profiled
# on a grid as in the spatial model. geoR is not used.

# Matern temporal correlation R_t(nu) and its derivative dR_t/dnu, using the same
# a = sqrt(2 kappa) |t - s| / nu convention as the spatial kernel.
.temporal_corr <- function(times, nu, kappa = 0.5, deriv = FALSE) {
  U <- abs(outer(times, times, "-"))
  cc <- sqrt(2 * kappa)
  a <- cc * U / nu
  e <- exp(-a)
  if (kappa == 0.5) {
    R <- e; dR <- a * e / nu
  } else if (kappa == 1.5) {
    R <- (1 + a) * e; dR <- a^2 * e / nu
  } else if (kappa == 2.5) {
    R <- (1 + a + a^2 / 3) * e; dR <- a^2 * (1 + a) * e / (3 * nu)
  } else stop("temporal kappa must be 0.5, 1.5 or 2.5")
  diag(R) <- 1; if (deriv) diag(dR) <- 0
  if (deriv) list(R = R, dR = dR) else R
}

# Kronecker-free per-sample joint log densities (vector over the nsim samples) for
# fixed spatial correlation Rs (at one phi). theta = (beta, log sigma2, log nu).
.st_num_loglik <- function(theta, D, Rsinv, ldetRs, S.sim, data_ll, N, T, p,
                            times, kappa_t, want_grad = FALSE) {
  beta <- theta[1:p]; s2 <- exp(theta[p + 1]); nu <- exp(theta[p + 2])
  tc <- .temporal_corr(times, nu, kappa_t, deriv = want_grad)
  Rt <- if (want_grad) tc$R else tc
  cRt <- chol(Rt); Rtinv <- chol2inv(cRt); ldetRt <- 2 * sum(log(diag(cRt)))
  mu <- as.numeric(D %*% beta)                      # length NT
  nsim <- nrow(S.sim)
  const <- N * T * log(s2) + N * ldetRt + T * ldetRs

  q <- numeric(nsim)
  Glist <- if (want_grad) vector("list", nsim) else NULL
  qnu <- if (want_grad) numeric(nsim) else NULL
  RtiDp <- if (want_grad) Rtinv %*% tc$dR %*% Rtinv else NULL
  for (b in seq_len(nsim)) {
    M <- matrix(S.sim[b, ] - mu, N, T)              # regions x times
    RsM <- Rsinv %*% M
    G <- RsM %*% Rtinv                              # Rs^{-1} M Rt^{-1}
    q[b] <- sum(G * M)
    if (want_grad) {
      Glist[[b]] <- G
      qnu[b] <- sum((RsM %*% RtiDp) * M)            # tr(Rs^{-1} M Rt^{-1} Rt' Rt^{-1} M')
    }
  }
  num <- data_ll - 0.5 * (const + q / s2)
  out <- list(num = num, q = q, s2 = s2, nu = nu)
  if (want_grad) { out$G <- Glist; out$qnu <- qnu; out$trRtiD <- sum(Rtinv * t(tc$dR)) }
  out
}

# Fit (beta, sigma2, nu) at one phi by maximising the MC log-likelihood.
.st_fit_one_phi <- function(theta0, D, Rsinv, ldetRs, S.sim, data_ll, Den,
                            N, T, p, times, kappa_t) {
  negMCL <- function(theta) {
    cp <- .st_num_loglik(theta, D, Rsinv, ldetRs, S.sim, data_ll, N, T, p, times, kappa_t)
    -log(mean(exp(cp$num - Den)))
  }
  negGrad <- function(theta) {
    cp <- .st_num_loglik(theta, D, Rsinv, ldetRs, S.sim, data_ll, N, T, p, times,
                         kappa_t, want_grad = TRUE)
    w <- exp(cp$num - Den); w <- w / sum(w)
    s2 <- cp$s2; nu <- cp$nu; nsim <- length(w)
    gb <- matrix(0, nsim, p)
    for (b in seq_len(nsim)) gb[b, ] <- as.numeric(crossprod(D, as.vector(cp$G[[b]]))) / s2
    gs <- -N * T / 2 + cp$q / (2 * s2)
    gnu <- nu * (-0.5 * N * cp$trRtiD + 0.5 * cp$qnu / s2)
    -colSums(w * cbind(gb, gs, gnu))
  }
  opt <- stats::nlminb(theta0, negMCL, negGrad)
  list(par = opt$par, value = -opt$objective, grad = negGrad)
}

#' Fit a spatio-temporal SDA-LGCP model (Kronecker-free)
#'
#' Separable space-time SDA-LGCP for aggregated counts observed over the same
#' \code{N} regions at \code{T} times. The spatial scale \code{phi} is profiled on
#' a grid; the temporal Matern range \code{nu} is estimated continuously. The
#' likelihood never forms the \eqn{(NT)\times(NT)} covariance.
#'
#' @param formula model formula (with optional \code{offset(log(pop))}).
#' @param data data frame of \code{N*T} rows ordered by time then region (rows
#'   \code{(t-1)*N + 1:N} are time \code{t}).
#' @param my_shp \code{sf} polygons for the \code{N} regions.
#' @param times numeric vector of length \code{T} of observation times.
#' @param delta candidate-point spacing.
#' @param phi spatial-scale grid (default from geometry).
#' @param kappa spatial Matern smoothness; \code{kappa_t} temporal smoothness.
#' @param method,weighted,pop_shp point-generation controls.
#' @param control.mcmc list from \code{\link{control_mcmc}}.
#' @param messages logical; print progress.
#' @return an object of class \code{c("SDALGCP2_ST","SDALGCP2")} with \code{phi_opt},
#'   \code{nu_opt}, coefficient table and covariance.
#' @export
SDALGCP2_ST <- function(formula, data, my_shp, times, delta, phi = NULL,
                        kappa = 0.5, kappa_t = 0.5, method = 3L, weighted = FALSE,
                        pop_shp = NULL, control.mcmc = NULL, messages = FALSE) {
  if (!inherits(my_shp, "sf")) my_shp <- sf::st_as_sf(my_shp)
  N <- nrow(my_shp); T <- length(times)
  mf <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mf))
  if (length(y) != N * T) stop("data must have N*T rows (N regions, T times).")
  D <- stats::model.matrix(attr(mf, "terms"), data); p <- ncol(D)
  m <- if (any(startsWith(names(mf), "offset"))) exp(stats::model.offset(mf)) else rep(1, N * T)

  if (is.null(phi)) {
    areas <- as.numeric(sf::st_area(my_shp)); bb <- sf::st_bbox(my_shp)
    phi <- seq(sqrt(min(areas)),
               min(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) / 10, length.out = 12)
  }
  if (is.null(control.mcmc)) control.mcmc <- control_mcmc(h = 1.65 / (N * T)^(1 / 6))

  pts <- sda_points(my_shp, delta, method = method, weighted = weighted, pop_shp = pop_shp)
  corr <- precompute_corr(pts, phi, kappa = kappa)
  Rarr <- corr$R

  # starting values from a Poisson GLM
  g <- stats::glm(formula, family = "poisson", data = data)
  beta0 <- stats::coef(g); s20 <- mean(stats::residuals(g)^2); nu0 <- stats::median(diff(sort(unique(times)))) * 2
  phi0 <- stats::median(phi); i0 <- which.min(abs(phi - phi0)); Rs0 <- Rarr[, , i0]
  Rt0 <- .temporal_corr(times, nu0, kappa_t)

  # one-time anchor sampling from the FULL covariance (Kronecker built once here)
  Sigma0 <- s20 * kronecker(Rt0, Rs0)
  mu0 <- as.numeric(D %*% beta0)
  if (messages) cat("Sampling latent field at anchor (NT =", N * T, ")...\n")
  S.sim <- laplace_sampling(mu0, Sigma0, y, m, control.mcmc)$samples
  data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)

  cRs0 <- chol(Rs0); Rs0inv <- chol2inv(cRs0); ldetRs0 <- 2 * sum(log(diag(cRs0)))
  theta0 <- c(beta0, log(s20), log(nu0))
  Den <- .st_num_loglik(theta0, D, Rs0inv, ldetRs0, S.sim, data_ll, N, T, p, times, kappa_t)$num

  if (messages) cat("Profiling over", length(phi), "phi values...\n")
  res <- vector("list", length(phi)); th <- theta0
  for (i in seq_along(phi)) {
    cRs <- chol(Rarr[, , i]); Rsinv <- chol2inv(cRs); ldetRs <- 2 * sum(log(diag(cRs)))
    res[[i]] <- .st_fit_one_phi(th, D, Rsinv, ldetRs, S.sim, data_ll, Den, N, T, p, times, kappa_t)
    th <- res[[i]]$par
  }
  vals <- vapply(res, `[[`, numeric(1), "value")
  best <- which.max(vals); est <- res[[best]]$par

  # covariance from a finite-difference Hessian of the analytic gradient
  gfun <- res[[best]]$grad
  eps <- 1e-4; k <- length(est); H <- matrix(0, k, k)
  g0 <- gfun(est)
  for (j in seq_len(k)) { e <- est; e[j] <- e[j] + eps; H[, j] <- (gfun(e) - g0) / eps }
  H <- (H + t(H)) / 2
  cov <- tryCatch(solve(H), error = function(e) matrix(NA, k, k))

  beta_opt <- est[1:p]; sigma2_opt <- exp(est[p + 1]); nu_opt <- exp(est[p + 2])
  pn <- c(colnames(D), "sigma^2", "nu")
  estimates <- stats::setNames(c(beta_opt, sigma2_opt, nu_opt), pn)
  dimnames(cov) <- list(pn, pn)

  out <- list(D = D, y = y, m = m, beta_opt = beta_opt, sigma2_opt = sigma2_opt,
              phi_opt = phi[best], nu_opt = nu_opt, estimates = estimates, cov = cov,
              llike_val_opt = vals[best], mu = mu0, kappa = kappa, kappa_t = kappa_t,
              all_para = data.frame(phi = phi, value = vals), S = S.sim,
              N = N, T = T, times = times, phi_method = "grid+direct(nu)", call = match.call())
  attr(out, "my_shp") <- my_shp
  class(out) <- c("SDALGCP2_ST", "SDALGCP2")
  out
}
