# Covariates measured on a DIFFERENT support from the outcome (e.g. pollution
# monitors at point locations). Two steps (see math/confounding-and-misalignment):
#   1. predict the covariate to the candidate points by Gaussian-process kriging,
#      giving a predictive mean mu_z and variance v_z at each point;
#   2. fit with the BERKSON-CORRECTED intensity-scale offset
#        b_i(beta) = log sum_k w_ik exp( z_tilde_ik' beta + 0.5 beta' V_ik beta ),
#      whose 0.5 beta' V_ik beta term propagates the prediction uncertainty (so beta
#      is not attenuated by plugging in the kriged mean). v_z -> 0 recovers the
#      raster case.

# --- Gaussian-process kriging (exponential covariance + nugget) -----------------
#
# The same machinery handles point and areal covariate data; only the data-data
# correlation R_BB(phi) and the data-target cross-correlation R_tB(phi) differ
# (plain exponential for points; aggregated over the covariate polygons for areal).

# Cross Euclidean distances between rows of A (m x 2) and B (n x 2).
.cross_dist <- function(A, B) {
  sqrt(outer(A[, 1], B[, 1], "-")^2 + outer(A[, 2], B[, 2], "-")^2)
}

# ML fit of a constant-mean GP given a builder RBB_fun(phi) for the L x L data
# correlation. Profiles out the mean and variance; optimises log phi and log nugget.
.fit_gp_generic <- function(z, RBB_fun, dmax) {
  z <- as.numeric(z); L <- length(z); one <- rep(1, L)
  negll <- function(par) {
    phi <- exp(par[1]); eta <- exp(par[2])
    C0 <- RBB_fun(phi); diag(C0) <- diag(C0) + eta
    ch <- tryCatch(chol(C0), error = function(e) NULL); if (is.null(ch)) return(1e10)
    Ci <- chol2inv(ch); ldet <- 2 * sum(log(diag(ch)))
    Ci1 <- Ci %*% one; mu <- sum(Ci1 * z) / sum(Ci1); r <- z - mu
    s2 <- as.numeric(crossprod(r, Ci %*% r)) / L
    0.5 * (L * log(s2) + ldet)
  }
  opt <- stats::optim(c(log(dmax / 3), log(0.1)), negll, method = "BFGS",
                      control = list(maxit = 200))
  phi <- exp(opt$par[1]); eta <- exp(opt$par[2])
  C0 <- RBB_fun(phi); diag(C0) <- diag(C0) + eta; ch <- chol(C0); Ci <- chol2inv(ch)
  Ci1 <- Ci %*% one; mu <- sum(Ci1 * z) / sum(Ci1); r <- z - mu
  s2 <- as.numeric(crossprod(r, Ci %*% r)) / L
  list(z = z, phi = phi, eta = eta, sigma2 = s2, mu = mu, Ci = Ci,
       Cir = as.numeric(Ci %*% r), Ci1 = as.numeric(Ci1), sumCi1 = sum(Ci1))
}

# Ordinary-kriging predictive mean/variance given the T x L cross-correlation RtB.
.krige_generic <- function(gp, RtB) {
  mu <- gp$mu + as.numeric(RtB %*% gp$Cir)
  quad <- rowSums((RtB %*% gp$Ci) * RtB)
  m1 <- 1 - as.numeric(RtB %*% gp$Ci1)
  v <- gp$sigma2 * pmax(0, 1 - quad + m1^2 / gp$sumCi1)
  list(mean = mu, var = v)
}

# Point covariate support: plain exponential correlation.
.gp_point <- function(s, z) {
  s <- as.matrix(s); D <- as.matrix(stats::dist(s))
  gp <- .fit_gp_generic(z, function(phi) exp(-D / phi), max(D)); gp$s <- s; gp
}
.krige_point <- function(gp, t)
  .krige_generic(gp, exp(-.cross_dist(as.matrix(t), gp$s) / gp$phi))

# Areal covariate support: the covariate is an average over polygons. Data-data and
# data-target covariances are the AGGREGATED exponential over candidate points
# inside the covariate polygons -- reusing the outcome-model kernels.
.gp_areal <- function(cov_sf, zcol, kappa = 0.5) {
  cpts <- sda_points(cov_sf, .auto_delta(cov_sf, 16), method = 3L)
  coords <- lapply(cpts, function(p) as.matrix(p$xy)[, 1:2, drop = FALSE])
  bb <- sf::st_bbox(cov_sf)
  dmax <- sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
  RBB <- function(phi) corr_aggregate_cpp(coords, list(), c(phi), kappa, FALSE, 0L)[, , 1]
  gp <- .fit_gp_generic(cov_sf[[zcol]], RBB, as.numeric(dmax))
  gp$coords <- coords; gp$kappa <- kappa; gp
}
.krige_areal <- function(gp, t)
  .krige_generic(gp, cross_cov_cpp(as.matrix(t), gp$coords, list(), gp$phi, gp$kappa, FALSE, 0L))

# --- Berkson-corrected tilting ---------------------------------------------------

# Like .tilt(), but each candidate point carries a covariate predictive variance
# Vlist[[i]] (n_i x q, columns matching Zlist; intercept column = 0). The exponent
# gains 0.5 beta' V beta and the effective design gains V beta.
.tilt_berkson <- function(Zlist, Vlist, wlist, beta) {
  N <- length(Zlist); q <- length(beta)
  Dc <- matrix(0, N, q); bvec <- numeric(N); clist <- vector("list", N)
  for (i in seq_len(N)) {
    Z <- Zlist[[i]]; V <- Vlist[[i]]; w <- wlist[[i]]
    # exponent: z'beta + 0.5 sum_j v_kj beta_j^2 ; gradient: z + (v_kj beta_j)_j
    berk <- 0.5 * as.numeric(V %*% (beta^2))
    lin  <- as.numeric(Z %*% beta) + berk
    mx <- max(lin); ew <- w * exp(lin - mx); sden <- sum(ew)
    bvec[i] <- mx + log(sden)
    ci <- ew / sden; clist[[i]] <- ci
    Veff <- V * rep(beta, each = nrow(V))           # column j scaled by beta_j  (= V beta)
    Dc[i, ] <- as.numeric(crossprod(Z + Veff, ci))  # sum_k c_ik (z_k + V_k beta)
  }
  colnames(Dc) <- colnames(Zlist[[1]])
  list(Dc = Dc, b = bvec, c = clist)
}

#' Fit an SDA-LGCP with covariates measured on a different support
#'
#' Covariates observed on a \emph{different support} from the outcome (e.g. air-
#' quality monitors at point locations) are kriged to the candidate points and
#' enter the model on the intensity scale with a Berkson correction that propagates
#' the prediction uncertainty (see \code{math/confounding-and-misalignment.pdf}).
#'
#' @param formula model formula; the covariate names appear on the right-hand side.
#' @param data \code{sf} polygons holding the response and offset (one row/region).
#' @param delta candidate-point spacing.
#' @param covariates a named list; each element is an \code{sf} carrying a column
#'   of the same name -- the covariate's observed values on its own support, either
#'   \strong{points} (e.g. monitors; plain kriging) or \strong{polygons} (areal
#'   averages on a different partition; aggregated areal kriging).
#' @param phi spatial-scale grid for the outcome model (default from geometry).
#' @param method,weighted,pop_shp point-generation controls.
#' @param berkson logical; include the Berkson uncertainty correction (default
#'   \code{TRUE}). \code{FALSE} gives the naive kriged-mean plug-in.
#' @param control.mcmc list from \code{\link{control_mcmc}}.
#' @param max_iter,tol outer Gauss-Newton controls.
#' @param messages logical; print progress.
#' @return an object of class \code{"SDALGCP2"} with \code{misaligned = TRUE}.
#' @seealso \code{\link{SDALGCP2_raster}}, \code{\link{sdalgcp}}
#' @export
SDALGCP2_misaligned <- function(formula, data, delta, covariates, phi = NULL,
                                method = 3L, weighted = FALSE, pop_shp = NULL,
                                berkson = TRUE, control.mcmc = NULL,
                                max_iter = 10L, tol = 1e-3, messages = FALSE) {
  if (!inherits(data, "sf")) data <- sf::st_as_sf(data)
  rhs <- attr(stats::terms(formula), "term.labels")
  rhs <- rhs[!grepl("^offset\\(", rhs)]
  miss <- setdiff(rhs, names(covariates))
  if (length(miss)) stop("no covariate data supplied for: ", paste(miss, collapse = ", "))
  for (nm in rhs) if (!nm %in% names(data)) data[[nm]] <- 0
  mf <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mf))
  m <- if (any(startsWith(names(mf), "offset"))) exp(stats::model.offset(mf)) else rep(1, length(y))

  if (is.null(phi)) {
    areas <- as.numeric(sf::st_area(data)); bb <- sf::st_bbox(data)
    phi <- seq(sqrt(min(areas)),
               min(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) / 10, length.out = 12)
  }

  pts <- sda_points(data, delta, method = method, weighted = weighted, pop_shp = pop_shp)
  wlist <- lapply(pts, function(p) as.numeric(p$weight))
  npt <- vapply(pts, function(p) nrow(as.matrix(p$xy)), 0L)
  allxy <- do.call(rbind, lapply(pts, function(p) as.matrix(p$xy)[, 1:2, drop = FALSE]))
  idx <- rep(seq_along(pts), npt)

  # krige each covariate from its own support to every candidate point
  q <- length(rhs) + 1L
  MU <- matrix(0, nrow(allxy), q); VAR <- matrix(0, nrow(allxy), q)
  MU[, 1] <- 1                                       # intercept column (no variance)
  for (j in seq_along(rhs)) {
    cv <- covariates[[rhs[j]]]
    if (!inherits(cv, "sf")) stop("covariate '", rhs[j], "' must be an sf object.")
    gtype <- as.character(sf::st_geometry_type(cv, by_geometry = FALSE))[1]
    areal <- grepl("POLYGON", gtype)
    if (messages) cat(sprintf("  kriging covariate '%s' from %d %s...\n",
                              rhs[j], nrow(cv), if (areal) "polygons" else "points"))
    if (areal) {
      gp <- .gp_areal(cv, rhs[j], kappa = 0.5); kr <- .krige_areal(gp, allxy)
    } else {
      gp <- .gp_point(sf::st_coordinates(cv)[, 1:2, drop = FALSE], cv[[rhs[j]]])
      kr <- .krige_point(gp, allxy)
    }
    MU[, j + 1] <- kr$mean; VAR[, j + 1] <- if (berkson) kr$var else 0
  }
  Zlist <- lapply(seq_along(pts), function(i) {
    Z <- MU[idx == i, , drop = FALSE]; colnames(Z) <- c("(Intercept)", rhs); Z })
  Vlist <- lapply(seq_along(pts), function(i) VAR[idx == i, , drop = FALSE])

  corr <- precompute_corr(pts, phi)
  beta <- stats::setNames(rep(0, q), c("(Intercept)", rhs))
  beta[1] <- log(mean(y / m))
  fit <- NULL; iter <- 0L
  repeat {
    iter <- iter + 1L
    tl <- .tilt_berkson(Zlist, Vlist, wlist, beta)
    off <- log(m) + tl$b - as.numeric(tl$Dc %*% beta)
    df <- as.data.frame(tl$Dc); names(df) <- paste0("Z", seq_len(q))
    df$.y <- y; df$.off <- off
    form <- stats::as.formula(paste0(".y ~ ",
              paste(paste0("Z", seq_len(q)), collapse = " + "), " - 1 + offset(.off)"))
    fit <- mcml_fit(form, df, corr, control.mcmc = control.mcmc, messages = messages)
    beta_new <- stats::setNames(fit$beta_opt, c("(Intercept)", rhs))
    db <- max(abs(beta_new - beta)); beta <- beta_new
    if (messages) cat(sprintf("  misaligned iter %d: max|dbeta| = %.4g\n", iter, db))
    if (db < tol || iter >= max_iter) break
  }
  names(fit$beta_opt) <- c("(Intercept)", rhs)
  pn <- c("(Intercept)", rhs, "sigma^2"); dimnames(fit$cov) <- list(pn, pn)
  fit$estimates <- stats::setNames(c(beta, fit$sigma2_opt), pn)
  fit$misaligned <- TRUE; fit$berkson <- berkson; fit$n_iter <- iter
  attr(fit, "my_shp") <- data; fit$call <- match.call()
  fit
}
