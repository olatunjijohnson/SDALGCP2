# Raster and misaligned (kriged) covariates for the spatio-temporal SDA-LGCP.
#
# As in the spatial case, the covariates are a time-invariant spatial surface that
# enters the intensity at the candidate-point level and is aggregated on the
# intensity (exp) scale via a log-sum-exp offset b_i(beta) -- not by averaging the
# predictor over each region. We linearise with the intensity-tilted effective
# design (the c_ik-weighted covariate) and fit the resulting linear model with the
# Kronecker-free spatio-temporal fitter SDALGCP2_ST(), iterating the tilting to
# convergence (Gauss-Newton). Because the covariate surface does not change over
# time, the per-region tilt is computed once and replicated across the T times;
# only the population offset varies by region and time.

# Gauss-Newton fit of an intensity-tilted covariate ST model. Zlist/Vlist/wlist are
# the per-region point covariates (Vlist NULL for rasters, the Berkson predictive
# variance for misaligned covariates); qn are the coefficient names. The inner
# fitter is SDALGCP2_ST() on the linearised (effective-design + offset) model.
.fit_st_tilted <- function(formula, data, my_shp, times, Zlist, Vlist, wlist, qn,
                           delta, phi, kappa, kappa_t, method, weighted, pop_shp,
                           control.mcmc, reanchor, max_iter, tol, messages) {
  N <- nrow(my_shp); T <- length(times); q <- length(qn)
  mf <- stats::model.frame(formula, data)
  y  <- as.numeric(stats::model.response(mf))
  if (length(y) != N * T) stop("data must have N*T rows (N regions, T times).")
  m  <- if (any(startsWith(names(mf), "offset"))) exp(stats::model.offset(mf)) else rep(1, N * T)
  # rows are time-blocked, region-fastest: row (t-1)*N + i is region i at time t
  region <- rep(seq_len(N), times = T)

  tilt <- function(beta) if (is.null(Vlist)) .tilt(Zlist, wlist, beta)
                         else .tilt_berkson(Zlist, Vlist, wlist, beta)

  beta <- stats::setNames(rep(0, q), qn); beta[1] <- log(mean(y / m))
  fit <- NULL; iter <- 0L
  repeat {
    iter <- iter + 1L
    tl <- tilt(beta)                                       # Dc (N x q), b (N)
    Dc_full <- tl$Dc[region, , drop = FALSE]               # (N*T) x q
    off <- log(m) + tl$b[region] - as.numeric(Dc_full %*% beta)
    df <- as.data.frame(Dc_full); names(df) <- paste0("Z", seq_len(q))
    df$.y <- y; df$.off <- off
    form <- stats::as.formula(paste0(".y ~ ",
              paste(paste0("Z", seq_len(q)), collapse = " + "), " - 1 + offset(.off)"))
    fit <- SDALGCP2_ST(form, df, my_shp, times = times, delta = delta, phi = phi,
                       kappa = kappa, kappa_t = kappa_t, method = method, weighted = weighted,
                       pop_shp = pop_shp, control.mcmc = control.mcmc, reanchor = reanchor,
                       messages = messages)
    beta_new <- stats::setNames(as.numeric(fit$beta_opt)[seq_len(q)], qn)
    db <- max(abs(beta_new - beta))
    if (messages) cat(sprintf("  ST-covariate iter %d: max|dbeta| = %.4g\n", iter, db))
    beta <- beta_new
    if (db < tol || iter >= max_iter) break
  }

  # relabel coefficients with the real covariate names and finalise
  names(fit$beta_opt) <- qn
  pn <- c(qn, "sigma^2", "nu")
  fit$estimates <- stats::setNames(c(fit$beta_opt, fit$sigma2_opt, fit$nu_opt), pn)
  dimnames(fit$cov) <- list(pn, pn)
  fit$n_iter <- iter
  fit
}

# Restricted spatial regression (RSR) for the spatio-temporal model -- addressing
# space-time confounding by the same device as the spatial case (Reich, Hodges &
# Zadnik 2006; Hughes & Haran 2013), now with the separable space-time covariance.
# With K an orthonormal basis of the orthogonal complement of the fixed-effect
# design D (over all N*T region-times),
#
#   Y ~ Poisson(m exp(eta)),  eta = D beta + K alpha,
#   alpha ~ N(0, sigma2 K'(Rt(nu) %x% Rs(phi))K),
#
# so the space-time random effect cannot reproduce anything in col(D) and beta is
# identified by the data, not the spatial-temporal structure. We integrate alpha
# out by a Laplace approximation (analytic mode), profile over the phi grid and
# estimate the temporal range nu. Reduces to the spatial restricted fit at T = 1.
# This forms the full (N*T) covariance, so it is O((N*T)^3) -- best for modest N*T.
.fit_restricted_st <- function(formula, data, my_shp, times, delta, phi = NULL,
                               kappa = 0.5, kappa_t = 0.5, method = 3L,
                               weighted = FALSE, pop_shp = NULL, messages = FALSE) {
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

  pts  <- sda_points(my_shp, delta, method = method, weighted = weighted, pop_shp = pop_shp)
  Rarr <- precompute_corr(pts, phi, kappa = kappa)$R          # N x N x nphi
  K <- .null_basis(D); r <- ncol(K)                           # (N*T) x (N*T - p)

  g <- stats::glm(formula, family = "poisson", data = data)
  beta0 <- stats::coef(g); s20 <- max(mean(stats::residuals(g)^2), 0.05)
  nu0 <- stats::median(diff(sort(unique(times)))); if (!is.finite(nu0) || nu0 <= 0) nu0 <- 1
  par0 <- c(beta0, log(s20), log(nu0))

  # K'(Rt(nu) %x% Rs)K with a small ridge (K'(.)K is near-singular at large phi).
  build_Sa <- function(Rs, nu) {
    Rt <- .temporal_corr(times, nu, kappa_t)
    Sa <- crossprod(K, kronecker(Rt, Rs) %*% K)
    Sa + diag(1e-6 * mean(diag(Sa)), nrow(Sa))
  }
  # Laplace-marginal log-likelihood at (beta, log sigma2, log nu) for a fixed Rs.
  marglik <- function(par, Rs) {
    beta <- par[1:p]; s2 <- exp(par[p + 1]); nu <- exp(par[p + 2])
    Sa <- build_Sa(Rs, nu); ch <- tryCatch(chol(Sa), error = function(e) NULL)
    if (is.null(ch)) return(-1e10)
    SaInv <- chol2inv(ch); ldetSa <- 2 * sum(log(diag(ch)))
    am <- .alpha_mode(y, m, as.numeric(D %*% beta), K, SaInv / s2)
    a <- am$alpha
    data_ll  <- sum(y * am$eta - m * exp(am$eta))
    prior_ll <- -0.5 * (r * log(s2) + ldetSa + as.numeric(crossprod(a, SaInv %*% a)) / s2)
    ldetH <- as.numeric(determinant(am$negH, logarithm = TRUE)$modulus)
    data_ll + prior_ll - 0.5 * ldetH
  }

  if (messages) cat("RSR (space-time): profiling", length(phi), "phi values...\n")
  res <- lapply(seq_along(phi), function(i) {
    Rs <- Rarr[, , i]
    opt <- stats::optim(par0, function(par) -marglik(par, Rs), method = "BFGS",
                        control = list(maxit = 200))
    list(par = opt$par, value = -opt$value)
  })
  vals <- vapply(res, `[[`, numeric(1), "value")
  best <- which.max(vals); est <- res[[best]]$par
  beta_opt <- est[1:p]; sigma2_opt <- exp(est[p + 1]); nu_opt <- exp(est[p + 2])
  Rs_b <- Rarr[, , best]

  nll <- function(par) -marglik(par, Rs_b)
  H <- .num_hessian(nll, est)
  cov <- tryCatch(solve(H), error = function(e) matrix(NA, p + 2, p + 2))
  J <- diag(c(rep(1, p), sigma2_opt, nu_opt))            # delta method (log -> natural)
  cov <- J %*% cov %*% t(J)
  pn <- c(colnames(D), "sigma^2", "nu"); dimnames(cov) <- list(pn, pn)

  Dbeta <- as.numeric(D %*% beta_opt)
  am <- .alpha_mode(y, m, Dbeta, K, chol2inv(chol(build_Sa(Rs_b, nu_opt))) / sigma2_opt)

  out <- list(D = D, y = y, m = m, beta_opt = beta_opt, sigma2_opt = sigma2_opt,
              nu_opt = nu_opt, phi_opt = phi[best],
              estimates = stats::setNames(c(beta_opt, sigma2_opt, nu_opt), pn),
              cov = cov, llike_val_opt = vals[best], mu = Dbeta,
              all_para = data.frame(phi = phi, value = vals),
              N = N, T = T, times = times, kappa = kappa, kappa_t = kappa_t,
              confounding = "restricted", K = K, alpha_hat = am$alpha, eta_hat = am$eta,
              S.coord = pts, control.mcmc = control_mcmc(h = 1.65 / (N * T)^(1 / 6)),
              phi_method = "grid+laplace(nu)", call = match.call())
  attr(out, "my_shp") <- my_shp
  class(out) <- c("SDALGCP2_ST", "SDALGCP2")
  out
}
