test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("GP kriging recovers a smooth field at held-out points", {
  set.seed(1)
  f <- function(xy) 1.5 * sin(xy[, 1] / 4) + cos(xy[, 2] / 5)
  s <- cbind(runif(60, 0, 20), runif(60, 0, 20)); z <- f(s)
  gp <- SDALGCP2:::.gp_point(s, z)
  t <- cbind(runif(30, 2, 18), runif(30, 2, 18))
  kr <- SDALGCP2:::.krige_point(gp, t)
  expect_lt(sqrt(mean((kr$mean - f(t))^2)), 0.5)     # predictions close to truth
  expect_true(all(kr$var >= 0))                      # non-negative kriging variance
})

test_that("Berkson tilt reduces to ordinary tilt when the variance is zero", {
  set.seed(2)
  Z <- list(cbind(`(Intercept)` = 1, z = rnorm(5)), cbind(`(Intercept)` = 1, z = rnorm(4)))
  V0 <- lapply(Z, function(z) matrix(0, nrow(z), ncol(z)))
  w  <- list(rep(0.2, 5), rep(0.25, 4))
  b  <- c(0.1, 0.7)
  a <- SDALGCP2:::.tilt_berkson(Z, V0, w, b)
  o <- SDALGCP2:::.tilt(Z, w, b)
  expect_equal(a$b, o$b, tolerance = 1e-12)
  expect_equal(a$Dc, o$Dc, tolerance = 1e-12)
})

test_that("misaligned-covariate fit recovers beta from point-observed covariate", {
  skip_on_cran()
  suppressMessages(library(sf))
  set.seed(123)
  zf <- function(xy) 1.5 * sin(xy[, 1] / 4) + 1.2 * cos(xy[, 2] / 5) + 0.5 * xy[, 1] / 20
  regions <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(7, 7)))
  N <- nrow(regions)
  pts <- sda_points(regions, delta = 1.2, method = 3); w <- lapply(pts, function(p) p$weight)
  bt <- c(-6, 0.8)
  b_true <- sapply(seq_len(N), function(i) {
    Z <- cbind(1, zf(as.matrix(pts[[i]]$xy))); log(sum(w[[i]] * exp(as.numeric(Z %*% bt)))) })
  regions$pop <- round(runif(N, 1000, 5000)); regions$cases <- rpois(N, regions$pop * exp(b_true))
  mon <- st_as_sf(data.frame(x = runif(40, 0, 20), y = runif(40, 0, 20)), coords = c("x", "y"))
  mon$z <- zf(st_coordinates(mon))
  ctrl <- control_mcmc(n.sim = 3000, burnin = 800, thin = 5, h = 1.65 / N^(1/6))

  fit <- SDALGCP2_misaligned(cases ~ z + offset(log(pop)), regions, delta = 1.2,
                             covariates = list(z = mon), control.mcmc = ctrl, max_iter = 6)
  expect_true(isTRUE(fit$misaligned))
  expect_lt(abs(fit$beta_opt["z"] - 0.8), 0.3)       # recovered from monitor data alone
  expect_true(is.finite(sqrt(fit$cov["z", "z"])))
})

test_that("areal misaligned covariate (different polygons) recovers beta", {
  skip_on_cran()
  suppressMessages(library(sf))
  set.seed(42)
  zf <- function(xy) 1.5 * sin(xy[, 1] / 4) + 1.2 * cos(xy[, 2] / 5) + 0.5 * xy[, 1] / 20
  regions <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(7, 7)))
  N <- nrow(regions)
  pts <- sda_points(regions, delta = 1.2, method = 3); w <- lapply(pts, function(p) p$weight)
  b_true <- sapply(seq_len(N), function(i) {
    Z <- cbind(1, zf(as.matrix(pts[[i]]$xy))); log(sum(w[[i]] * exp(as.numeric(Z %*% c(-6, 0.8))))) })
  regions$pop <- round(runif(N, 1000, 5000)); regions$cases <- rpois(N, regions$pop * exp(b_true))
  # covariate observed as averages over a coarser, different 4x4 partition
  covpoly <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(4, 4)))
  cpts <- sda_points(covpoly, delta = 1.2, method = 3)
  covpoly$z <- sapply(cpts, function(p) mean(zf(as.matrix(p$xy))))
  ctrl <- control_mcmc(n.sim = 3000, burnin = 800, thin = 5, h = 1.65 / N^(1/6))
  fit <- SDALGCP2_misaligned(cases ~ z + offset(log(pop)), regions, delta = 1.2,
                             covariates = list(z = covpoly), control.mcmc = ctrl, max_iter = 6)
  expect_true(isTRUE(fit$misaligned))
  expect_lt(abs(fit$beta_opt["z"] - 0.8), 0.35)
})
