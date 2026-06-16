# Spatially continuous (raster) covariates for SDA-LGCP (task #5).
#
# Covariates z(x) supplied as rasters enter the LGCP intensity at the POINT level,
# inside the exponential, and are aggregated on the intensity scale -- NOT by
# averaging the predictor over each polygon (which is biased under the nonlinear
# log link). The aggregated Poisson mean is
#
#   E[Y_i | T_i] = m_i * exp( b_i(beta) + T_i ),   T_i ~ N(0, sigma2 R^w(phi)),
#   b_i(beta) = log sum_k w_ik exp( z(x_ik)' beta )          (log-sum-exp offset)
#
# i.e. covariates enter through the log-sum-exp of the within-region covariate
# surface, b_i(beta) != (zbar_i)'beta. We fit by a Gauss-Newton fixed point: the
# local linearisation of b_i has gradient db_i/dbeta = sum_k c_ik z_ik, the
# intensity-tilted ("effective") covariate, with tilted weights
# c_ik = w_ik exp(z_ik'beta) / sum_l w_il exp(z_il'beta). Iterating a linear MCML
# fit with this effective design and a matched offset converges to the MLE of the
# nonlinear model, and reuses mcml_fit() unchanged.

# Extract raster covariate values at every candidate point of every region.
# Returns a list of n_i x q matrices (with an intercept column prepended).
.point_covariates <- function(points, rasters, crs) {
  if (!inherits(rasters, "SpatRaster")) rasters <- terra::rast(rasters)
  lapply(points, function(p) {
    xy <- as.matrix(p$xy)[, 1:2, drop = FALSE]
    v <- terra::extract(rasters, xy)
    v <- as.matrix(v[, !(names(v) %in% "ID"), drop = FALSE])
    storage.mode(v) <- "double"
    cbind(`(Intercept)` = 1, v)
  })
}

# One intensity-scale tilting step: returns the effective design D^c (N x q),
# the log-sum-exp offsets b (length N) and the tilted weights per region.
.tilt <- function(Zlist, wlist, beta) {
  N <- length(Zlist); q <- length(beta)
  Dc <- matrix(0, N, q); bvec <- numeric(N); clist <- vector("list", N)
  for (i in seq_len(N)) {
    Z <- Zlist[[i]]; w <- wlist[[i]]
    lin <- as.numeric(Z %*% beta)
    mx <- max(lin)                       # numerical stability for log-sum-exp
    ew <- w * exp(lin - mx)
    sden <- sum(ew)
    bvec[i] <- mx + log(sden)            # b_i(beta)
    ci <- ew / sden                      # tilted weights, sum_k = 1
    clist[[i]] <- ci
    Dc[i, ] <- as.numeric(crossprod(Z, ci))   # sum_k c_ik z_ik  (db_i/dbeta)
  }
  colnames(Dc) <- colnames(Zlist[[1]])
  list(Dc = Dc, b = bvec, c = clist)
}

#' Fit an SDA-LGCP with spatially continuous (raster) covariates
#'
#' Covariates supplied as rasters enter the model at the candidate-point level and
#' are aggregated on the intensity (exp) scale via a log-sum-exp offset
#' \eqn{b_i(\beta)=\log\sum_k w_{ik}\exp(z(x_{ik})^\top\beta)} -- the statistically
#' correct alternative to averaging the predictor over each polygon. Estimation is
#' a Gauss-Newton fixed point that reuses \code{\link{mcml_fit}} with the
#' intensity-tilted effective design.
#'
#' @param formula model formula; right-hand-side names must match raster layer
#'   names. The response and an \code{offset(log(pop))} come from \code{data}.
#' @param data data frame with the response and offset (one row per region).
#' @param my_shp \code{sf} polygons.
#' @param delta candidate-point spacing.
#' @param rasters a \code{terra::SpatRaster} (or object coercible by
#'   \code{terra::rast}) whose layers are the spatially varying covariates named
#'   in \code{formula}.
#' @param phi spatial-scale grid (default chosen from the geometry).
#' @param method,weighted,pop_shp point-generation controls (see \code{\link{sda_points}}).
#' @param control.mcmc list from \code{\link{control_mcmc}}.
#' @param max_iter,tol outer Gauss-Newton iteration controls.
#' @param messages logical; print progress.
#' @return an object of class \code{"SDALGCP2"} (as \code{\link{mcml_fit}}) with
#'   extra fields \code{raster = TRUE} and \code{n_iter}.
#' @seealso \code{\link{SDALGCP2}}, \code{\link{mcml_fit}}
#' @export
SDALGCP2_raster <- function(formula, data, my_shp, delta, rasters, phi = NULL,
                            method = 3L, weighted = FALSE, pop_shp = NULL,
                            control.mcmc = NULL, max_iter = 10L, tol = 1e-3,
                            messages = FALSE) {
  if (!inherits(my_shp, "sf")) my_shp <- sf::st_as_sf(my_shp)
  if (!inherits(rasters, "SpatRaster")) rasters <- terra::rast(rasters)

  # response and offset (population) from the formula/data
  mf <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mf))
  m <- if (any(startsWith(names(mf), "offset"))) exp(stats::model.offset(mf)) else rep(1, length(y))
  rhs <- attr(stats::terms(formula), "term.labels")
  rhs <- rhs[!grepl("^offset\\(", rhs)]
  miss <- setdiff(rhs, names(rasters))
  if (length(miss)) stop("raster layers missing for covariate(s): ", paste(miss, collapse = ", "))

  if (is.null(phi)) {
    areas <- as.numeric(sf::st_area(my_shp)); bb <- sf::st_bbox(my_shp)
    phi <- seq(sqrt(min(areas)),
               min(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) / 10, length.out = 15)
  }

  pts <- sda_points(my_shp, delta, method = method, weighted = weighted, pop_shp = pop_shp)
  Zlist <- .point_covariates(pts, rasters[[rhs]], sf::st_crs(my_shp))
  wlist <- lapply(pts, function(p) as.numeric(p$weight))
  q <- ncol(Zlist[[1]]); qn <- colnames(Zlist[[1]])

  corr <- precompute_corr(pts, phi)              # R^w(phi): does not depend on beta

  beta <- stats::setNames(rep(0, q), qn)
  beta[1] <- log(mean(y / m))                    # intercept start
  fit <- NULL; iter <- 0L
  repeat {
    iter <- iter + 1L
    tl <- .tilt(Zlist, wlist, beta)
    off <- log(m) + tl$b - as.numeric(tl$Dc %*% beta)   # linearisation offset
    df <- as.data.frame(tl$Dc); names(df) <- paste0("Z", seq_len(q))
    df$.y <- y; df$.off <- off
    form <- stats::as.formula(paste0(".y ~ ",
              paste(paste0("Z", seq_len(q)), collapse = " + "), " - 1 + offset(.off)"))
    fit <- mcml_fit(form, df, corr, control.mcmc = control.mcmc, messages = messages)
    beta_new <- stats::setNames(fit$beta_opt, qn)
    delta_b <- max(abs(beta_new - beta))
    if (messages) cat(sprintf("  raster iter %d: max|dbeta| = %.4g\n", iter, delta_b))
    beta <- beta_new
    if (delta_b < tol || iter >= max_iter) break
  }

  # relabel coefficients with the real covariate names and finalise
  names(fit$beta_opt) <- qn
  pn <- c(qn, "sigma^2")
  dimnames(fit$cov) <- list(pn, pn)
  fit$estimates <- stats::setNames(c(beta, fit$sigma2_opt), pn)
  fit$raster <- TRUE; fit$n_iter <- iter; fit$call <- match.call()
  fit
}
