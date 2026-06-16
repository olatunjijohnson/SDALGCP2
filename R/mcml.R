# Vectorised Monte Carlo maximum likelihood for the spatial SDA-LGCP (bottleneck B3).
# Replaces SDALGCP's Aggregated_poisson_log_MCML()/SDALGCPParaEst() apply/Reduce
# loops with Cholesky-based matrix operations. Statistically identical: same
# importance-sampling MCML objective, gradient and Hessian.

# Per-sample log joint density components for a fixed phi (R fixed via Rinv/ldetR).
# Returns a list of reusable pieces; `data_ll` (sum y*S - m*exp(S)) is independent
# of the parameters and is supplied precomputed.
.mc_num_loglik <- function(par, D, Rinv, ldetR, S.sim, data_ll, n, p) {
  beta <- par[1:p]
  sigma2 <- exp(par[p + 1])
  mu <- as.numeric(D %*% beta)
  diff <- S.sim - rep(mu, each = nrow(S.sim))          # nsim x n
  G <- diff %*% Rinv                                   # nsim x n
  Q <- rowSums(G * diff)                               # quadratic forms
  num <- data_ll - 0.5 * (n * log(sigma2) + ldetR + Q / sigma2)
  list(num = num, diff = diff, G = G, Q = Q, mu = mu, sigma2 = sigma2)
}

# Fit (beta, sigma2) for one fixed phi by maximising the MC log-likelihood.
.mcml_one_phi <- function(par0_opt, D, Rinv, ldetR, S.sim, data_ll,
                          Denominator, n, p, DtRinvD, DtRinv, messages = FALSE) {
  nsim <- nrow(S.sim)

  neg_mcl <- function(par) {
    cp <- .mc_num_loglik(par, D, Rinv, ldetR, S.sim, data_ll, n, p)
    -log(mean(exp(cp$num - Denominator)))
  }

  neg_grad <- function(par) {
    cp <- .mc_num_loglik(par, D, Rinv, ldetR, S.sim, data_ll, n, p)
    w <- exp(cp$num - Denominator); w <- w / sum(w)
    GB <- (cp$G %*% D) / cp$sigma2                      # nsim x p: rows are grad.beta
    gs <- -n / 2 + 0.5 * cp$Q / cp$sigma2               # grad wrt log(sigma2)
    grad_beta <- colSums(w * GB)
    grad_s <- sum(w * gs)
    -c(grad_beta, grad_s)
  }

  neg_hess <- function(par) {
    cp <- .mc_num_loglik(par, D, Rinv, ldetR, S.sim, data_ll, n, p)
    w <- exp(cp$num - Denominator); w <- w / sum(w)
    GB <- (cp$G %*% D) / cp$sigma2                      # nsim x p
    gs <- -n / 2 + 0.5 * cp$Q / cp$sigma2               # nsim
    grad_i <- cbind(GB, gs)                             # nsim x (p+1)
    gbar <- colSums(w * grad_i)

    GG <- crossprod(grad_i * w, grad_i)                 # sum_i w_i grad grad'
    # mean Hessian of the per-sample joint density (sum_i w_i Hess L_i):
    meanH <- matrix(0, p + 1, p + 1)
    meanH[1:p, 1:p] <- -DtRinvD / cp$sigma2
    hbs <- -colSums(w * GB)
    meanH[1:p, p + 1] <- hbs
    meanH[p + 1, 1:p] <- hbs
    meanH[p + 1, p + 1] <- -0.5 * sum(w * cp$Q) / cp$sigma2

    H <- meanH + GG - tcrossprod(gbar)
    -H
  }

  opt <- stats::nlminb(par0_opt, neg_mcl, neg_grad, neg_hess,
                       control = list(trace = as.integer(messages)))
  H <- neg_hess(opt$par)
  list(par = opt$par, value = -opt$objective, cov = solve(H))
}

#' Monte Carlo maximum likelihood estimation for the spatial SDA-LGCP
#'
#' Vectorised, Cholesky-based re-implementation of \code{SDALGCP::SDALGCPParaEst()}.
#' Simulates the latent field at an anchor, then profiles the importance-sampling
#' MCML objective over the supplied \code{phi} grid.
#'
#' @param formula model formula, optionally with an \code{offset()} term.
#' @param data data frame holding the model variables.
#' @param corr list with \code{R} (N x N x n_phi correlation array) and \code{phi},
#'   e.g. from \code{\link{precompute_corr}}.
#' @param par0 optional starting values \code{c(beta, sigma2, phi)}; if \code{NULL}
#'   they are derived from a Poisson GLM.
#' @param control.mcmc list from \code{\link{control_mcmc}} (defaults if \code{NULL}).
#' @param messages logical; print optimiser progress.
#' @return an object of class \code{"SDALGCP2"} (estimates, covariance, profile,
#'   latent samples and metadata).
#' @export
mcml_fit <- function(formula, data, corr, par0 = NULL, control.mcmc = NULL,
                     phi_method = c("grid", "direct"), nugget = FALSE,
                     reanchor = 0L, reanchor_tol = 1e-2, messages = FALSE) {
  phi_method <- match.arg(phi_method)
  if (nugget && phi_method != "direct")
    stop("nugget = TRUE requires phi_method = 'direct'.")
  mf <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mf))
  D <- stats::model.matrix(attr(mf, "terms"), data)
  n <- length(y); p <- ncol(D)
  m <- if (any(startsWith(names(mf), "offset"))) exp(stats::model.offset(mf)) else rep(1, n)

  phi <- as.numeric(corr$phi)
  R <- corr$R
  if (is.null(control.mcmc)) control.mcmc <- control_mcmc(h = 1.65 / n^(1 / 6))

  if (is.null(par0)) {
    g <- stats::glm(formula, family = "poisson", data = data)
    par0 <- c(stats::coef(g), mean(stats::residuals(g)^2), stats::median(phi))
  }
  if (any(par0[-(1:p)] <= 0)) stop("Covariance parameters in 'par0' must be positive.")

  # A single MCML pass: draw the latent field at the supplied anchor and fit.
  # Wrapped so re-anchoring can re-run it with the previous optimum as the anchor.
  .pass <- function(par0, nu_anchor = 0.1) {
  beta0 <- par0[1:p]; sigma2_0 <- par0[p + 1]; phi0 <- par0[p + 2]

  # Anchor correlation = grid entry nearest phi0 (plus the nugget, so the latent
  # field is drawn from the full covariance the importance sampler reweights to).
  i0 <- which.min(abs(phi - phi0))
  R0 <- R[, , i0]
  if (nugget) diag(R0) <- diag(R0) + nu_anchor
  Sigma0 <- sigma2_0 * R0
  mu0 <- as.numeric(D %*% beta0)

  S.sim <- laplace_sampling(mu0, Sigma0, y, m, control.mcmc)$samples

  # Pieces independent of phi/par.
  data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)

  ## ---- continuous-phi ("direct") path: optimise (beta, sigma2, phi) jointly ----
  if (phi_method == "direct") {
    if (!is.null(corr$kappa) && corr$kappa != 0.5)
      stop("phi_method = 'direct' currently supports the exponential kernel (kappa = 0.5).")
    pts <- attr(R, "S_coord")
    weighted <- isTRUE(attr(R, "weighted"))
    coords <- lapply(pts, function(z) as.matrix(z$xy)[, 1:2, drop = FALSE])
    wts <- if (weighted) lapply(pts, function(z) as.numeric(z$weight)) else list()
    if (nugget) {
      df <- .mcml_direct_nugget_fit(y, D, m, coords, wts, weighted, S.sim, data_ll,
                                    par0_opt = c(beta0, log(sigma2_0)), phi0 = phi0,
                                    nu0 = nu_anchor, n = n, p = p, messages = messages)
      pnames <- c(colnames(D), "sigma^2", "phi", "nu")
      phi_opt <- df$estimate[p + 2]; nu_opt <- df$estimate[p + 3]
      Cm <- corr_and_grad_cpp(coords, wts, phi_opt, weighted, 0L)$R
      diag(Cm) <- diag(Cm) + nu_opt
    } else {
      df <- .mcml_direct_fit(y, D, m, coords, wts, weighted, S.sim, data_ll,
                             par0_opt = c(beta0, log(sigma2_0)), phi0 = phi0,
                             n = n, p = p, messages = messages)
      pnames <- c(colnames(D), "sigma^2", "phi")
      phi_opt <- df$estimate[p + 2]; nu_opt <- NULL
      Cm <- corr_and_grad_cpp(coords, wts, phi_opt, weighted, 0L)$R
    }
    beta_opt <- df$estimate[1:p]; sigma2_opt <- df$estimate[p + 1]
    out <- list(
      D = D, y = y, m = m,
      beta_opt = beta_opt, sigma2_opt = sigma2_opt, phi_opt = phi_opt, nu_opt = nu_opt,
      estimates = stats::setNames(df$estimate, pnames),
      cov = df$cov, Sigma_mat_opt = sigma2_opt * Cm,
      llike_val_opt = df$value, mu = as.numeric(D %*% beta_opt),
      all_para = data.frame(phi = phi_opt, value = df$value),
      all_cov = list(df$cov), par0 = par0, control.mcmc = control.mcmc,
      S = S.sim, S.coord = pts,
      kappa = if (!is.null(corr$kappa)) corr$kappa else 0.5,
      phi_method = "direct", nugget = nugget, call = NULL
    )
    attr(out, "weighted") <- weighted
    attr(out, "my_shp") <- attr(R, "my_shp")
    attr(out, "S_coord") <- pts
    attr(out, "prematrix") <- corr
    class(out) <- "SDALGCP2"
    return(out)
  }

  # Denominator (importance anchor): per-sample log joint at (beta0, sigma2_0, R0).
  ch0 <- chol(R0); R0inv <- chol2inv(ch0); ldetR0 <- 2 * sum(log(diag(ch0)))
  Denominator <- .mc_num_loglik(c(beta0, log(sigma2_0)), D, R0inv, ldetR0,
                                S.sim, data_ll, n, p)$num

  par0_opt <- c(beta0, log(sigma2_0))
  np <- length(phi)
  res <- vector("list", np)
  for (i in seq_len(np)) {
    ch <- chol(R[, , i]); Rinv <- chol2inv(ch); ldetR <- 2 * sum(log(diag(ch)))
    DtRinv <- crossprod(D, Rinv); DtRinvD <- DtRinv %*% D
    res[[i]] <- .mcml_one_phi(par0_opt, D, Rinv, ldetR, S.sim, data_ll,
                              Denominator, n, p, DtRinvD, DtRinv, messages)
    par0_opt <- res[[i]]$par                          # warm start next phi
  }

  vals <- vapply(res, `[[`, numeric(1), "value")
  best <- which.max(vals)
  est <- res[[best]]$par
  beta_opt <- est[1:p]; sigma2_opt <- exp(est[p + 1])
  pnames <- c(colnames(D), "sigma^2")
  for (i in seq_along(res)) dimnames(res[[i]]$cov) <- list(pnames, pnames)

  all_para <- data.frame(phi = phi, value = vals,
                         t(vapply(res, function(z) c(z$par[1:p], exp(z$par[p + 1])),
                                  numeric(p + 1))))
  colnames(all_para) <- c("phi", "value", colnames(D), "sigma2")

  out <- list(
    D = D, y = y, m = m,
    beta_opt = beta_opt, sigma2_opt = sigma2_opt, phi_opt = phi[best],
    estimates = stats::setNames(c(beta_opt, sigma2_opt), pnames),
    cov = res[[best]]$cov,
    Sigma_mat_opt = sigma2_opt * R[, , best],
    llike_val_opt = vals[best],
    mu = as.numeric(D %*% beta_opt),
    all_para = all_para, all_cov = lapply(res, `[[`, "cov"),
    par0 = par0, control.mcmc = control.mcmc, S = S.sim,
    S.coord = attr(R, "S_coord"), kappa = if (!is.null(corr$kappa)) corr$kappa else 0.5,
    phi_method = "grid", call = NULL
  )
  attr(out, "weighted") <- attr(R, "weighted")
  attr(out, "my_shp") <- attr(R, "my_shp")
  attr(out, "S_coord") <- attr(R, "S_coord")
  attr(out, "prematrix") <- corr
  class(out) <- "SDALGCP2"
  out
  }  # end .pass

  # Re-anchoring loop: refit using the previous optimum as the new anchor so the
  # importance weights stay near-uniform (raises the MC effective sample size).
  cur <- par0; nu_cur <- 0.1; fit <- NULL; done <- 0L
  for (k in 0:reanchor) {
    fit_new <- .pass(cur, nu_anchor = nu_cur)
    conv <- FALSE
    if (k > 0) {
      prev <- c(fit$beta_opt, fit$sigma2_opt, fit$phi_opt)
      neww <- c(fit_new$beta_opt, fit_new$sigma2_opt, fit_new$phi_opt)
      conv <- max(abs(neww - prev) / pmax(abs(prev), 1e-6)) < reanchor_tol
    }
    fit <- fit_new; done <- k
    cur <- c(fit$beta_opt, fit$sigma2_opt, fit$phi_opt)
    if (nugget && !is.null(fit$nu_opt)) nu_cur <- max(fit$nu_opt, 1e-3)
    if (messages && reanchor > 0) cat(sprintf("  re-anchor pass %d done\n", k))
    if (conv) break
  }
  fit$n_reanchor <- done
  fit$call <- match.call()
  fit
}
