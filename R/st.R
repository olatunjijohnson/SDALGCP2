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
#' @param kappa spatial Matern smoothness.
#' @param kappa_t temporal Matern smoothness.
#' @param method,weighted,pop_shp point-generation controls.
#' @param control.mcmc list from \code{\link{control_mcmc}}.
#' @param reanchor number of re-anchoring passes (re-simulate the latent field at
#'   the current optimum and refit); improves the variance-parameter estimates.
#' @param messages logical; print progress.
#' @return an object of class \code{c("SDALGCP2_ST","SDALGCP2")} with \code{phi_opt},
#'   \code{nu_opt}, coefficient table and covariance.
#' @seealso \code{\link{sdalgcp}} (friendly wrapper), \code{\link{predict.SDALGCP2_ST}}
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' shp <- sdalgcp_data
#' ## build a 3-time panel (data frame, N*T rows ordered by time then region)
#' times <- 1:3
#' dat <- do.call(rbind, lapply(times, function(t) {
#'   d <- sf::st_drop_geometry(shp); d$time <- t
#'   d$cases <- rpois(nrow(d), d$pop * exp(-6 + 0.6 * d$x1 + 0.1 * (t - 2)))
#'   d
#' }))
#' fit <- SDALGCP2_ST(cases ~ x1 + offset(log(pop)), dat, shp, times = times,
#'                    delta = 1.5,
#'                    control.mcmc = control_mcmc(n.sim = 2000, burnin = 500, thin = 5))
#' fit$phi_opt; fit$nu_opt
#' }
#' @export
SDALGCP2_ST <- function(formula, data, my_shp, times, delta, phi = NULL,
                        kappa = 0.5, kappa_t = 0.5, method = 3L, weighted = FALSE,
                        pop_shp = NULL, control.mcmc = NULL, reanchor = 0L,
                        messages = FALSE) {
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
  beta0 <- stats::coef(g); s20 <- mean(stats::residuals(g)^2)
  nu0 <- stats::median(diff(sort(unique(times)))) * 2
  mu0 <- as.numeric(D %*% beta0)

  # A single pass: draw the latent field at the anchor, profile over phi.
  .pass <- function(beta_a, s2_a, nu_a, phi_a) {
    i0 <- which.min(abs(phi - phi_a)); Rs0 <- Rarr[, , i0]
    Rt0 <- .temporal_corr(times, nu_a, kappa_t)
    Sigma0 <- s2_a * kronecker(Rt0, Rs0)            # full Kronecker built ONCE per pass
    if (messages) cat("Sampling latent field at anchor (NT =", N * T, ")...\n")
    S.sim <- laplace_sampling(as.numeric(D %*% beta_a), Sigma0, y, m, control.mcmc)$samples
    data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)
    cRs0 <- chol(Rs0); Rs0inv <- chol2inv(cRs0); ldetRs0 <- 2 * sum(log(diag(cRs0)))
    theta0 <- c(beta_a, log(s2_a), log(nu_a))
    Den <- .st_num_loglik(theta0, D, Rs0inv, ldetRs0, S.sim, data_ll, N, T, p, times, kappa_t)$num

    res <- vector("list", length(phi)); th <- theta0
    for (i in seq_along(phi)) {
      cRs <- chol(Rarr[, , i]); Rsinv <- chol2inv(cRs); ldetRs <- 2 * sum(log(diag(cRs)))
      res[[i]] <- .st_fit_one_phi(th, D, Rsinv, ldetRs, S.sim, data_ll, Den, N, T, p, times, kappa_t)
      th <- res[[i]]$par
    }
    vals <- vapply(res, `[[`, numeric(1), "value")
    best <- which.max(vals); est <- res[[best]]$par

    gfun <- res[[best]]$grad
    eps <- 1e-4; k <- length(est); H <- matrix(0, k, k); g0 <- gfun(est)
    for (j in seq_len(k)) { e <- est; e[j] <- e[j] + eps; H[, j] <- (gfun(e) - g0) / eps }
    H <- (H + t(H)) / 2
    cov <- tryCatch(solve(H), error = function(e) matrix(NA, k, k))
    list(beta_opt = est[1:p], sigma2_opt = exp(est[p + 1]), nu_opt = exp(est[p + 2]),
         phi_opt = phi[best], cov = cov, value = vals[best],
         all_para = data.frame(phi = phi, value = vals), S = S.sim)
  }

  # Re-anchoring loop (re-simulate the latent field at the current optimum).
  cur <- list(beta_opt = beta0, sigma2_opt = s20, nu_opt = nu0, phi_opt = stats::median(phi))
  fit_p <- NULL; done <- 0L
  for (kk in 0:reanchor) {
    fit_p <- .pass(cur$beta_opt, cur$sigma2_opt, cur$nu_opt, cur$phi_opt)
    conv <- kk > 0 && max(abs(c(fit_p$beta_opt, fit_p$sigma2_opt, fit_p$nu_opt) -
                              c(cur$beta_opt, cur$sigma2_opt, cur$nu_opt)) /
                          pmax(abs(c(cur$beta_opt, cur$sigma2_opt, cur$nu_opt)), 1e-6)) < 1e-2
    cur <- fit_p; done <- kk
    if (messages && reanchor > 0) cat(sprintf("  ST re-anchor pass %d done\n", kk))
    if (conv) break
  }

  pn <- c(colnames(D), "sigma^2", "nu")
  estimates <- stats::setNames(c(fit_p$beta_opt, fit_p$sigma2_opt, fit_p$nu_opt), pn)
  dimnames(fit_p$cov) <- list(pn, pn)
  out <- list(D = D, y = y, m = m, beta_opt = fit_p$beta_opt, sigma2_opt = fit_p$sigma2_opt,
              phi_opt = fit_p$phi_opt, nu_opt = fit_p$nu_opt, estimates = estimates,
              cov = fit_p$cov, llike_val_opt = fit_p$value, mu = mu0,
              kappa = kappa, kappa_t = kappa_t, all_para = fit_p$all_para, S = fit_p$S,
              S.coord = pts, control.mcmc = control.mcmc,
              N = N, T = T, times = times, n_reanchor = done,
              phi_method = "grid+direct(nu)", call = match.call())
  attr(out, "my_shp") <- my_shp
  class(out) <- c("SDALGCP2_ST", "SDALGCP2")
  out
}

#' Discrete (region x time) prediction for a spatio-temporal fit
#'
#' Draws the latent field at the fitted optimum and returns posterior mean and SD
#' of the incidence relative risk \eqn{\exp(\mu+S)} and covariate-adjusted relative
#' risk \eqn{\exp(S)} for every region and time.
#'
#' @param object an \code{"SDALGCP2_ST"} fit.
#' @param control.mcmc optional MCMC controls (defaults to the fitting ones).
#' @param ... unused.
#' @return a long \code{\link[sf]{sf}} of class
#'   \code{c("SDALGCP2_ST_pred", "sf", "data.frame")} with one row per region and
#'   time (ordered region-fastest within each time block) and columns
#'   \code{region}, \code{time}, \code{relative_risk}, \code{relative_risk_se}
#'   (\eqn{\exp(\mu+S)}), \code{adjusted_rr} and \code{adjusted_rr_se}
#'   (\eqn{\exp(S)}) -- the same column names as the spatial
#'   \code{\link{predict.SDALGCP2}}. The posterior draws are kept in object
#'   attributes (for exceedance); map a time slice with
#'   \code{\link{plot.SDALGCP2_ST_pred}}.
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' ## stack the spatial example into a 3-time panel with a mild temporal trend
#' times <- 1:3
#' panel <- do.call(rbind, lapply(times, function(t) {
#'   d <- sdalgcp_data; d$time <- t
#'   d$cases <- rpois(nrow(d), d$pop * exp(-6 + 0.6 * d$x1 + 0.1 * (t - 2)))
#'   d
#' }))
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = panel, time = "time",
#'                control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
#'                                          reanchor = 0))
#' pr <- predict(fit)        # a long sf: region x time
#' head(pr)
#' plot(pr, time = 2)        # map the relative risk at time 2
#' }
#' @method predict SDALGCP2_ST
#' @export
predict.SDALGCP2_ST <- function(object, control.mcmc = NULL, ...) {
  N <- object$N; T <- object$T
  pts <- object$S.coord
  if (is.null(pts)) stop("Prediction needs the candidate points stored in the fit.")
  Rs <- precompute_corr(pts, object$phi_opt, kappa = object$kappa)$R[, , 1]
  Rt <- .temporal_corr(object$times, object$nu_opt, object$kappa_t)
  Sigma <- object$sigma2_opt * kronecker(Rt, Rs)
  if (is.null(control.mcmc)) control.mcmc <- object$control.mcmc
  if (is.null(control.mcmc)) control.mcmc <- control_mcmc(h = 1.65 / (N * T)^(1 / 6))
  S.sim <- laplace_sampling(object$mu, Sigma, object$y, object$m, control.mcmc)$samples
  RR <- exp(S.sim); ARR <- exp(sweep(S.sim, 2, object$mu, "-"))

  # Long sf: one row per (region, time), region-fastest within each time block,
  # matching the column ordering of S.sim (vec of the N x T region-by-time field).
  shp  <- attr(object, "my_shp")
  geom <- sf::st_geometry(shp)[rep(seq_len(N), T)]
  out  <- sf::st_sf(region = rep(seq_len(N), T),
                    time   = rep(object$times, each = N),
                    relative_risk    = colMeans(RR), relative_risk_se = .colSDs(RR),
                    adjusted_rr      = colMeans(ARR), adjusted_rr_se   = .colSDs(ARR),
                    geometry = geom)
  attr(out, "pred_type")  <- "spatio-temporal"
  attr(out, "pred_draws") <- list(eta = S.sim, mu = object$mu)
  attr(out, "N") <- N; attr(out, "T") <- T; attr(out, "times") <- object$times
  class(out) <- c("SDALGCP2_ST_pred", class(out))
  out
}

#' Map a spatio-temporal prediction for one time
#'
#' Maps a chosen quantity (\code{"relative_risk"}, \code{"adjusted_rr"},
#' \code{"relative_risk_se"}, \code{"adjusted_rr_se"} or \code{"exceedance"}) for a
#' selected time slice of a spatio-temporal prediction.
#'
#' @param x an \code{"SDALGCP2_ST_pred"} object from \code{predict()} on an
#'   \code{"SDALGCP2_ST"} fit.
#' @param time the time to map (one of the fitted \code{times}); defaults to the
#'   first. Use \code{NULL} to facet all times.
#' @param what one of \code{"relative_risk"}, \code{"adjusted_rr"},
#'   \code{"relative_risk_se"}, \code{"adjusted_rr_se"}, \code{"exceedance"}.
#' @param threshold threshold for \code{what = "exceedance"}.
#' @param which for exceedance: \code{"adjusted_rr"} (default) or
#'   \code{"relative_risk"}.
#' @param ... unused.
#' @return a \code{ggplot} object.
#' @examples
#' \donttest{
#' data(sdalgcp_data)
#' times <- 1:3
#' panel <- do.call(rbind, lapply(times, function(t) {
#'   d <- sdalgcp_data; d$time <- t
#'   d$cases <- rpois(nrow(d), d$pop * exp(-6 + 0.6 * d$x1 + 0.1 * (t - 2)))
#'   d
#' }))
#' fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = panel, time = "time",
#'                control = sdalgcp_control(n_sim = 2000, burnin = 500, thin = 5,
#'                                          reanchor = 0))
#' pr <- predict(fit)
#' plot(pr, time = 2)                 # one time slice
#' plot(pr, time = NULL)              # facet all times
#' plot(pr, what = "exceedance", threshold = 1.2, time = 3)
#' }
#' @method plot SDALGCP2_ST_pred
#' @export
plot.SDALGCP2_ST_pred <- function(x, time = attr(x, "times")[1],
                                  what = c("relative_risk", "adjusted_rr",
                                           "relative_risk_se", "adjusted_rr_se", "exceedance"),
                                  threshold = 1, which = c("adjusted_rr", "relative_risk"), ...) {
  what <- match.arg(what); which <- match.arg(which)
  N <- attr(x, "N"); times_all <- attr(x, "times"); draws <- attr(x, "pred_draws")
  times <- if (is.null(time)) times_all else time
  cols <- match(times, times_all)
  if (anyNA(cols)) stop("'time' must be among the fitted times: ", paste(times_all, collapse = ", "))

  excd <- function(tcol) {
    d <- if (which == "adjusted_rr") exp(sweep(draws$eta, 2, draws$mu, "-")) else exp(draws$eta)
    idx <- ((tcol - 1) * N + 1):(tcol * N)
    apply(d[, idx, drop = FALSE], 2, function(v) mean(v > threshold))
  }
  getvals <- function(tcol) {
    rows <- ((tcol - 1) * N + 1):(tcol * N)
    switch(what,
      relative_risk = x$relative_risk[rows], adjusted_rr = x$adjusted_rr[rows],
      relative_risk_se = x$relative_risk_se[rows], adjusted_rr_se = x$adjusted_rr_se[rows],
      exceedance = excd(tcol))
  }

  maps <- do.call(rbind, lapply(cols, function(tc) {
    rows <- ((tc - 1) * N + 1):(tc * N)
    g <- x[rows, ]; g$fillvalue <- getvals(tc); g$.time <- times_all[tc]
    g[, c("fillvalue", ".time")] }))
  lab <- if (what == "exceedance") sprintf("P(%s > %g)", .which_label(which), threshold) else .var_label(what)
  mid <- if (what %in% c("relative_risk", "adjusted_rr")) 1 else NULL
  p <- ggplot2::ggplot(maps) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data$fillvalue), color = "grey70", linewidth = 0.1)
  if (length(cols) > 1) p <- p + ggplot2::facet_wrap(~ .time)
  sc <- if (what == "exceedance")
    ggplot2::scale_fill_viridis_c(name = lab, limits = c(0, 1), option = "magma")
  else if (!is.null(mid))
    ggplot2::scale_fill_gradient2(name = lab, midpoint = mid, low = "#2166AC",
                                  mid = "grey95", high = "#B2182B")
  else ggplot2::scale_fill_viridis_c(name = lab)
  p + sc + ggplot2::theme_minimal() +
    ggplot2::labs(title = if (length(cols) == 1) paste(what, "-", times) else what,
                  x = NULL, y = NULL)
}
