test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("mc_diagnostics returns a sane ESS and MC standard error", {
  suppressMessages(library(sf))
  set.seed(5)
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 16, ymax = 16))), n = c(6, 6)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 1.4, method = 3)
  phi_grid <- seq(1.5, 5, length.out = 6)
  corr <- precompute_corr(pts, phi_grid)
  Sig <- 0.4 * corr$R[, , 3]
  x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
  pop <- round(runif(N, 500, 3000))
  y <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.4) +
                          as.numeric(t(chol(Sig)) %*% rnorm(N))))
  dat <- data.frame(y = y, x1 = x1, pop = pop)
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / N^(1/6))
  fit <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.4,
                  phi = phi_grid, method = 3, control.mcmc = ctrl)

  d <- suppressWarnings(mc_diagnostics(fit))
  expect_true(d$ESS >= 1 && d$ESS <= d$B)
  expect_true(is.finite(d$se_loglik) && d$se_loglik >= 0)
  expect_equal(d$ESS_frac, d$ESS / d$B)
})

test_that("re-anchoring raises the importance-sampling ESS", {
  suppressMessages(library(sf))
  set.seed(8)
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 16, ymax = 16))), n = c(6, 6)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 1.4, method = 3)
  phi_grid <- seq(1, 6, length.out = 8)
  corr <- precompute_corr(pts, phi_grid)
  Sig <- 0.6 * corr$R[, , 6]                       # anchor at median will be far off
  x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
  pop <- round(runif(N, 500, 3000))
  y <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
                          as.numeric(t(chol(Sig)) %*% rnorm(N))))
  dat <- data.frame(y = y, x1 = x1, pop = pop)
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / N^(1/6))

  set.seed(1); f0 <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.4,
                              phi = phi_grid, method = 3, control.mcmc = ctrl, reanchor = 0)
  set.seed(1); f3 <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.4,
                              phi = phi_grid, method = 3, control.mcmc = ctrl, reanchor = 3)
  e0 <- suppressWarnings(mc_diagnostics(f0))$ESS_frac
  e3 <- suppressWarnings(mc_diagnostics(f3))$ESS_frac
  expect_true(f3$n_reanchor >= 1)
  expect_gt(e3, e0)             # re-anchoring improves the effective sample size
  expect_gt(e3, 0.5)           # and brings it to a healthy level
})
