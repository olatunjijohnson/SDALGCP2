# Continuous-phi ("direct") Monte Carlo maximum likelihood (task #6).
# Instead of profiling phi on a grid, optimise (beta, log sigma2, log phi) jointly.
# The scale phi sits inside the aggregated double integral R(phi); its first and
# second derivatives R'(phi), R''(phi) come from corr_and_grad_cpp and feed the
# analytic gradient/Hessian of the MC log-likelihood. Exponential kernel only.
#
# Notation (per importance sample b): u_b = S_b - D beta, q_b = u_b' Rinv u_b.
# Gradient of the per-sample joint log-density g_b w.r.t. xi = (beta, log s2, log phi):
#   d/dbeta   = D' Rinv u_b / s2
#   d/dlog s2 = -n/2 + q_b/(2 s2)
#   d/dlog phi= phi * ( -0.5 tr(Rinv R') + 0.5 (u_b' Rinv R' Rinv u_b)/s2 )
# MC log-lik grad = sum_b w_b grad g_b ; Hessian = sum_b w_b (Hess g_b + grad grad')
#                   - (sum_b w_b grad)(sum_b w_b grad)'.

# Build R(phi), derivatives and the theta-only matrices needed by grad/Hessian.
.phi_pieces <- function(phi, coords, wts, weighted, kappa = 0.5, need_hess = FALSE) {
  cg <- corr_and_grad_cpp(coords, wts, phi, kappa, weighted, 0L)
  R <- cg$R; Rp <- cg$dR
  ch <- chol(R); Rinv <- chol2inv(ch); ldetR <- 2 * sum(log(diag(ch)))
  A  <- Rinv %*% Rp                 # Rinv R'
  m2 <- A %*% Rinv                  # Rinv R' Rinv
  trA <- sum(diag(A))
  out <- list(R = R, Rinv = Rinv, ldetR = ldetR, A = A, m2 = m2, trA = trA)
  if (need_hess) {
    Rpp <- cg$d2R
    B2 <- Rinv %*% Rpp
    trB <- sum(diag(B2))
    trA2 <- sum(A * t(A))           # tr(A A)
    n2 <- Rinv %*% (2 * Rp %*% Rinv %*% Rp - Rpp) %*% Rinv
    out$n2 <- n2
    out$t2 <- -0.5 * (trB - trA2)
  }
  out
}

# Continuous-(phi, nu) direct fit WITH a relative nugget: covariance
# C(phi, nu) = R(phi) + nu I (so the full covariance is sigma2 (R + nu I)).
# Parameters xi = (beta, log sigma2, log phi, log nu); dC/dphi = R'(phi),
# dC/dnu = I, d2C/dnu2 = d2C/dphi dnu = 0. Analytic gradient + Hessian, gated by
# numerical differentiation (test-nugget.R).
.mcml_direct_nugget_fit <- function(y, D, m, coords, wts, weighted, S.sim, data_ll,
                                    par0_opt, phi0, nu0, n, p, kappa = 0.5,
                                    messages = FALSE) {
  B <- nrow(S.sim); s2idx <- p + 1; pidx <- p + 2; nidx <- p + 3

  pieces <- function(phi, nu) {
    cg <- corr_and_grad_cpp(coords, wts, phi, kappa, weighted, 0L)
    Rp <- cg$dR; Rpp <- cg$d2R
    C <- cg$R; diag(C) <- diag(C) + nu
    ch <- chol(C); Cinv <- chol2inv(ch); ldetC <- 2 * sum(log(diag(ch)))
    Aphi <- Cinv %*% Rp; m2phi <- Aphi %*% Cinv
    Cinv2 <- Cinv %*% Cinv
    n2phi <- Cinv %*% (2 * Rp %*% Cinv %*% Rp - Rpp) %*% Cinv
    list(C = C, Cinv = Cinv, ldetC = ldetC, Rp = Rp,
         Aphi = Aphi, m2phi = m2phi, trAphi = sum(diag(Aphi)),
         tr_m2phi = sum(diag(m2phi)),
         t2phi = -0.5 * (sum(Cinv * Rpp) - sum(Aphi * t(Aphi))),
         n2phi = n2phi, Cinv2 = Cinv2, trCinv = sum(diag(Cinv)),
         trCinv2 = sum(Cinv * Cinv), RpCinv = Rp %*% Cinv)
  }

  num_at <- function(beta, s2, pc) {
    mu <- as.numeric(D %*% beta)
    U <- S.sim - rep(mu, each = B)
    Q <- rowSums((U %*% pc$Cinv) * U)
    list(num = data_ll - 0.5 * (n * log(s2) + pc$ldetC + Q / s2), U = U, Q = Q)
  }

  pc0 <- pieces(phi0, nu0)
  Den <- num_at(par0_opt[1:p], exp(par0_opt[p + 1]), pc0)$num
  unpack <- function(xi) list(beta = xi[1:p], s2 = exp(xi[s2idx]),
                              phi = exp(xi[pidx]), nu = exp(xi[nidx]))

  negMCL <- function(xi) {
    th <- unpack(xi); pc <- pieces(th$phi, th$nu); cp <- num_at(th$beta, th$s2, pc)
    -log(mean(exp(cp$num - Den)))
  }

  grad_parts <- function(th, pc, cp) {
    s2 <- th$s2; phi <- th$phi; nu <- th$nu
    U <- cp$U; W <- U %*% pc$Cinv
    GB  <- (W %*% D) / s2
    gs  <- -n / 2 + cp$Q / (2 * s2)
    qm2phi <- rowSums((U %*% pc$m2phi) * U)
    qm2nu  <- rowSums(W * W)                       # u' Cinv^2 u
    gphi <- phi * (-0.5 * pc$trAphi + 0.5 * qm2phi / s2)
    gnu  <- nu  * (-0.5 * pc$trCinv + 0.5 * qm2nu  / s2)
    list(W = W, GB = GB, gs = gs, qm2phi = qm2phi, qm2nu = qm2nu,
         gphi = gphi, gnu = gnu)
  }

  negGrad <- function(xi) {
    th <- unpack(xi); pc <- pieces(th$phi, th$nu); cp <- num_at(th$beta, th$s2, pc)
    w <- exp(cp$num - Den); w <- w / sum(w)
    gpz <- grad_parts(th, pc, cp)
    -colSums(w * cbind(gpz$GB, gpz$gs, gpz$gphi, gpz$gnu))
  }

  negHess <- function(xi) {
    th <- unpack(xi); pc <- pieces(th$phi, th$nu); cp <- num_at(th$beta, th$s2, pc)
    w <- exp(cp$num - Den); w <- w / sum(w)
    s2 <- th$s2; phi <- th$phi; nu <- th$nu
    U <- cp$U; gpz <- grad_parts(th, pc, cp); W <- gpz$W
    grad_i <- cbind(gpz$GB, gpz$gs, gpz$gphi, gpz$gnu)
    gbar <- colSums(w * grad_i); GG <- crossprod(grad_i * w, grad_i)

    qn2phi <- rowSums((U %*% pc$n2phi) * U)
    qn2nu  <- 2 * rowSums((W %*% pc$Cinv) * W)     # u' (2 Cinv^3) u
    qphinu <- rowSums((W %*% pc$RpCinv) * W)       # w' R' Cinv w,  w = Cinv u
    DtCinvD <- crossprod(D, pc$Cinv) %*% D
    Bphi <- (U %*% pc$m2phi) %*% D                 # rows (D' m2phi u_b)'
    Bnu  <- (W %*% pc$Cinv) %*% D                  # rows (D' Cinv^2 u_b)'

    P <- p + 3; mH <- matrix(0, P, P)
    mH[1:p, 1:p] <- -DtCinvD / s2
    hbs <- -colSums(w * gpz$GB);             mH[1:p, s2idx] <- mH[s2idx, 1:p] <- hbs
    hbp <- -phi * colSums(w * Bphi) / s2;    mH[1:p, pidx]  <- mH[pidx, 1:p]  <- hbp
    hbn <- -nu  * colSums(w * Bnu)  / s2;    mH[1:p, nidx]  <- mH[nidx, 1:p]  <- hbn
    mH[s2idx, s2idx] <- -sum(w * cp$Q) / (2 * s2)
    mH[s2idx, pidx]  <- mH[pidx, s2idx] <- -phi * sum(w * gpz$qm2phi) / (2 * s2)
    mH[s2idx, nidx]  <- mH[nidx, s2idx] <- -nu  * sum(w * gpz$qm2nu)  / (2 * s2)
    mH[pidx, pidx] <- sum(w * (phi^2 * (pc$t2phi - 0.5 * qn2phi / s2) + gpz$gphi))
    mH[nidx, nidx] <- sum(w * (nu^2  * (0.5 * pc$trCinv2 - 0.5 * qn2nu / s2) + gpz$gnu))
    mH[pidx, nidx] <- mH[nidx, pidx] <-
      phi * nu * (0.5 * pc$tr_m2phi - sum(w * qphinu) / s2)
    -(mH + GG - tcrossprod(gbar))
  }

  xi0 <- c(par0_opt[1:p], par0_opt[p + 1], log(phi0), log(nu0))
  opt <- stats::nlminb(xi0, negMCL, negGrad, negHess, control = list(trace = as.integer(messages)))
  H <- negHess(opt$par); cov <- solve(H)
  est <- c(opt$par[1:p], exp(opt$par[s2idx]), exp(opt$par[pidx]), exp(opt$par[nidx]))
  pnames <- c(colnames(D), "sigma^2", "phi", "nu")
  dimnames(cov) <- list(pnames, pnames)
  list(estimate = est, value = -opt$objective, cov = cov,
       obj = negMCL, grad = negGrad, hess = negHess)
}

# One direct fit: returns estimate c(beta, sigma2, phi), value, covariance.
.mcml_direct_fit <- function(y, D, m, coords, wts, weighted, S.sim, data_ll,
                             par0_opt, phi0, n, p, kappa = 0.5, messages = FALSE) {
  B <- nrow(S.sim)

  num_at <- function(beta, s2, pc) {
    mu <- as.numeric(D %*% beta)
    diff <- S.sim - rep(mu, each = B)
    Q <- rowSums((diff %*% pc$Rinv) * diff)
    list(num = data_ll - 0.5 * (n * log(s2) + pc$ldetR + Q / s2),
         diff = diff, Q = Q, mu = mu)
  }

  # Anchor denominator at phi0 (prior part; data term cancels in the ratio).
  pc0 <- .phi_pieces(phi0, coords, wts, weighted, kappa)
  Den <- num_at(par0_opt[1:p], exp(par0_opt[p + 1]), pc0)$num

  unpack <- function(xi) list(beta = xi[1:p], s2 = exp(xi[p + 1]), phi = exp(xi[p + 2]))

  negMCL <- function(xi) {
    th <- unpack(xi); pc <- .phi_pieces(th$phi, coords, wts, weighted, kappa)
    cp <- num_at(th$beta, th$s2, pc)
    -log(mean(exp(cp$num - Den)))
  }

  negGrad <- function(xi) {
    th <- unpack(xi); pc <- .phi_pieces(th$phi, coords, wts, weighted, kappa)
    cp <- num_at(th$beta, th$s2, pc)
    w <- exp(cp$num - Den); w <- w / sum(w)
    G  <- cp$diff %*% pc$Rinv
    GB <- (G %*% D) / th$s2                       # B x p : d/dbeta
    gs <- -n / 2 + cp$Q / (2 * th$s2)             # d/dlog s2
    qm2 <- rowSums((cp$diff %*% pc$m2) * cp$diff)
    gp <- th$phi * (-0.5 * pc$trA + 0.5 * qm2 / th$s2)   # d/dlog phi
    -colSums(w * cbind(GB, gs, gp))
  }

  negHess <- function(xi) {
    th <- unpack(xi); pc <- .phi_pieces(th$phi, coords, wts, weighted, kappa, need_hess = TRUE)
    cp <- num_at(th$beta, th$s2, pc)
    w <- exp(cp$num - Den); w <- w / sum(w)
    s2 <- th$s2; phi <- th$phi
    G  <- cp$diff %*% pc$Rinv
    GB <- (G %*% D) / s2
    gs <- -n / 2 + cp$Q / (2 * s2)
    qm2 <- rowSums((cp$diff %*% pc$m2) * cp$diff)
    qn2 <- rowSums((cp$diff %*% pc$n2) * cp$diff)
    gp <- phi * (-0.5 * pc$trA + 0.5 * qm2 / s2)
    grad_i <- cbind(GB, gs, gp)                   # B x (p+2)
    gbar <- colSums(w * grad_i)
    GG <- crossprod(grad_i * w, grad_i)

    DtRinvD <- crossprod(D, pc$Rinv) %*% D
    Dm2 <- crossprod(D, pc$m2)                    # p x n  (D' m2)
    Bp <- cp$diff %*% t(Dm2)                      # B x p : (Dm2 u_b)'

    P <- p + 2
    meanH <- matrix(0, P, P)
    meanH[1:p, 1:p] <- -DtRinvD / s2
    hbs <- -colSums(w * GB)
    meanH[1:p, p + 1] <- meanH[p + 1, 1:p] <- hbs
    hbp <- -phi * colSums(w * Bp) / s2
    meanH[1:p, p + 2] <- meanH[p + 2, 1:p] <- hbp
    meanH[p + 1, p + 1] <- -sum(w * cp$Q) / (2 * s2)
    meanH[p + 1, p + 2] <- meanH[p + 2, p + 1] <- -phi * sum(w * qm2) / (2 * s2)
    meanH[p + 2, p + 2] <- sum(w * (phi^2 * (pc$t2 - 0.5 * qn2 / s2) + gp))

    -(meanH + GG - tcrossprod(gbar))
  }

  xi0 <- c(par0_opt[1:p], par0_opt[p + 1], log(phi0))
  opt <- stats::nlminb(xi0, negMCL, negGrad, negHess,
                       control = list(trace = as.integer(messages)))
  H <- negHess(opt$par)
  cov <- solve(H)
  est <- c(opt$par[1:p], exp(opt$par[p + 1]), exp(opt$par[p + 2]))
  pnames <- c(colnames(D), "sigma^2", "phi")
  dimnames(cov) <- list(pnames, pnames)
  list(estimate = est, value = -opt$objective, cov = cov,
       obj = negMCL, grad = negGrad, hess = negHess)   # closures exposed for validation
}
