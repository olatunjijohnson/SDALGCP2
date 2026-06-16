# Restricted spatial regression (RSR) for the SDA-LGCP, addressing spatial
# confounding (Reich, Hodges & Zadnik 2006; Hughes & Haran 2013). A spatially
# smooth covariate is collinear with the spatial random effect, which then absorbs
# the covariate's signal and biases beta. RSR constrains the random effect to the
# orthogonal complement of the fixed-effect design: with K an orthonormal basis of
# null(D'), the model is
#
#   Y_i ~ Poisson(m_i exp(eta_i)),   eta = D beta + K alpha,   alpha ~ N(0, sigma2 K'R(phi)K),
#
# so K alpha cannot reproduce anything in col(D) and beta is identified by the data,
# not by the spatial structure. We fit by a Laplace-approximate marginal likelihood
# (the latent alpha is integrated out analytically at its mode), profiling over phi.

# Orthonormal basis of the orthogonal complement of col(D)  (N x (N - p)).
.null_basis <- function(D) {
  qrd <- qr(D)
  Q <- qr.Q(qrd, complete = TRUE)
  Q[, (ncol(D) + 1):nrow(D), drop = FALSE]
}

# Mode of the restricted random effect alpha given (beta, sigma2, phi). Newton on
#   f(alpha) = sum(y eta - m e^eta) - 0.5 alpha' (Sa^{-1}/sigma2) alpha,  eta = Dbeta + K alpha.
.alpha_mode <- function(y, m, Dbeta, K, SaInv_s2, maxit = 100, tol = 1e-9) {
  alpha <- rep(0, ncol(K))
  for (it in seq_len(maxit)) {
    eta <- Dbeta + as.numeric(K %*% alpha)
    h <- m * exp(eta)
    grad <- as.numeric(crossprod(K, y - h)) - as.numeric(SaInv_s2 %*% alpha)
    negH <- crossprod(K * h, K) + SaInv_s2          # K' diag(h) K + Sa^{-1}/sigma2
    step <- solve(negH, grad)
    alpha <- alpha + step
    if (max(abs(step)) < tol) break
  }
  eta <- Dbeta + as.numeric(K %*% alpha)
  list(alpha = alpha, eta = eta, negH = crossprod(K * (m * exp(eta)), K) + SaInv_s2)
}

# Laplace-approximate marginal log-likelihood at (beta, log sigma2) for a fixed phi
# (with Sa = K'R(phi)K precomputed: SaInv, ldetSa).
.rsr_marglik <- function(par, y, m, D, K, SaInv, ldetSa, p, r) {
  beta <- par[1:p]; s2 <- exp(par[p + 1])
  Dbeta <- as.numeric(D %*% beta)
  am <- .alpha_mode(y, m, Dbeta, K, SaInv / s2)
  a <- am$alpha
  data_ll  <- sum(y * am$eta - m * exp(am$eta))
  prior_ll <- -0.5 * (r * log(s2) + ldetSa + as.numeric(crossprod(a, SaInv %*% a)) / s2)
  ldetH <- as.numeric(determinant(am$negH, logarithm = TRUE)$modulus)
  data_ll + prior_ll - 0.5 * ldetH                  # + constant (dropped)
}

# Restricted fit, profiling phi over the grid in `corr`. Returns an SDALGCP2 object.
.fit_restricted <- function(formula, data, corr, messages = FALSE) {
  mf <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mf))
  D <- stats::model.matrix(attr(mf, "terms"), data)
  n <- length(y); p <- ncol(D)
  m <- if (any(startsWith(names(mf), "offset"))) exp(stats::model.offset(mf)) else rep(1, n)
  phi <- as.numeric(corr$phi); R <- corr$R
  K <- .null_basis(D); r <- ncol(K)

  g <- stats::glm(formula, family = "poisson", data = data)
  par0 <- c(stats::coef(g), log(max(mean(stats::residuals(g)^2), 0.05)))

  fit_phi <- function(i) {
    Sa <- crossprod(K, R[, , i] %*% K)              # K' R(phi_i) K
    ch <- chol(Sa); SaInv <- chol2inv(ch); ldetSa <- 2 * sum(log(diag(ch)))
    opt <- stats::optim(par0, function(par) -.rsr_marglik(par, y, m, D, K, SaInv, ldetSa, p, r),
                        method = "BFGS", control = list(maxit = 200))
    list(par = opt$par, value = -opt$value, SaInv = SaInv, ldetSa = ldetSa)
  }

  if (messages) cat("RSR: profiling", length(phi), "phi values...\n")
  res <- lapply(seq_along(phi), fit_phi)
  vals <- vapply(res, `[[`, numeric(1), "value")
  best <- which.max(vals); rb <- res[[best]]
  est <- rb$par; beta_opt <- est[1:p]; sigma2_opt <- exp(est[p + 1])

  # Covariance of (beta, log sigma2) from the numerical Hessian of the profile.
  nll <- function(par) -.rsr_marglik(par, y, m, D, K, rb$SaInv, rb$ldetSa, p, r)
  H <- .num_hessian(nll, est)
  cov <- tryCatch(solve(H), error = function(e) matrix(NA, p + 1, p + 1))
  # delta-method: var(sigma2) = sigma2^2 var(log sigma2)
  J <- diag(c(rep(1, p), sigma2_opt)); cov <- J %*% cov %*% t(J)
  pn <- c(colnames(D), "sigma^2"); dimnames(cov) <- list(pn, pn)

  # latent mode at the optimum (for prediction)
  Dbeta <- as.numeric(D %*% beta_opt)
  am <- .alpha_mode(y, m, Dbeta, K, rb$SaInv / sigma2_opt)

  out <- list(D = D, y = y, m = m, beta_opt = beta_opt, sigma2_opt = sigma2_opt,
              phi_opt = phi[best], estimates = stats::setNames(c(beta_opt, sigma2_opt), pn),
              cov = cov, llike_val_opt = vals[best], mu = Dbeta,
              all_para = data.frame(phi = phi, value = vals),
              S = NULL, S.coord = attr(R, "S_coord"),
              kappa = if (!is.null(corr$kappa)) corr$kappa else 0.5,
              confounding = "restricted", K = K, alpha_hat = am$alpha,
              Sigma_mat_opt = sigma2_opt * R[, , best],
              eta_hat = am$eta, phi_method = "grid", call = match.call())
  attr(out, "weighted") <- attr(R, "weighted")
  attr(out, "my_shp") <- attr(R, "my_shp")
  attr(out, "S_coord") <- attr(R, "S_coord")
  attr(out, "prematrix") <- corr
  class(out) <- "SDALGCP2"
  out
}

# Simple finite-difference Hessian (central differences).
.num_hessian <- function(fn, x, eps = 1e-4) {
  k <- length(x); H <- matrix(0, k, k); f0 <- fn(x)
  for (i in seq_len(k)) for (j in i:k) {
    xi <- x; xi[i] <- xi[i] + eps; xi[j] <- xi[j] + eps; fpp <- fn(xi)
    xi <- x; xi[i] <- xi[i] + eps; xi[j] <- xi[j] - eps; fpm <- fn(xi)
    xi <- x; xi[i] <- xi[i] - eps; xi[j] <- xi[j] + eps; fmp <- fn(xi)
    xi <- x; xi[i] <- xi[i] - eps; xi[j] <- xi[j] - eps; fmm <- fn(xi)
    H[i, j] <- H[j, i] <- (fpp - fpm - fmp + fmm) / (4 * eps^2)
  }
  H
}
