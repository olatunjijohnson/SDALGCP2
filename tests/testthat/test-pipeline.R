test_that("testthat is operational", {
  expect_equal(1, 1)
})

make_sim <- function(seed = 7, ncell = 6) {
  suppressMessages(library(sf))
  set.seed(seed)
  bound <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 18, ymax = 18)))
  shp   <- st_sf(geometry = st_make_grid(bound, n = c(ncell, ncell)))
  N     <- nrow(shp)
  pts   <- sda_points(shp, delta = 1.5, method = 3)
  phi_grid <- seq(1, 5, length.out = 6)
  corr  <- precompute_corr(pts, phi_grid)
  Sig   <- 0.5 * corr$R[, , which.min(abs(phi_grid - 2.5))]
  x1    <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
  pop   <- round(runif(N, 500, 3000))
  y     <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
                              as.numeric(t(chol(Sig)) %*% rnorm(N))))
  list(shp = shp, dat = data.frame(y = y, x1 = x1, pop = pop),
       phi_grid = phi_grid, N = N)
}

test_that("sda_points returns weights summing to one per region", {
  suppressMessages(library(sf))
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 10))), n = c(4, 4)))
  pts <- sda_points(shp, delta = 1.5, method = 3)
  expect_length(pts, nrow(shp))
  expect_true(all(vapply(pts, function(p) abs(sum(p$weight) - 1) < 1e-8, logical(1))))
  expect_false(isTRUE(attr(pts, "weighted")))
})

test_that("SDALGCP2 end-to-end fit runs and recovers the slope", {
  sim <- make_sim()
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / sim$N^(1/6))
  fit <- SDALGCP2(y ~ x1 + offset(log(pop)), sim$dat, sim$shp, delta = 1.5,
                  phi = sim$phi_grid, method = 3, control.mcmc = ctrl)
  expect_s3_class(fit, "SDALGCP2")
  expect_lt(abs(fit$beta_opt[2] - 0.5), 0.3)
  expect_true(fit$sigma2_opt > 0)
})

test_that("prediction (discrete + continuous, both samplers) produces sane output", {
  sim <- make_sim()
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / sim$N^(1/6))
  fit <- SDALGCP2(y ~ x1 + offset(log(pop)), sim$dat, sim$shp, delta = 1.5,
                  phi = sim$phi_grid, method = 3, control.mcmc = ctrl)

  pd <- predict(fit, type = "discrete", sampler = "mcmc", control.mcmc = ctrl)
  expect_length(pd$RR_mean, sim$N)
  expect_true(all(pd$RR_se >= 0))

  pc <- predict(fit, type = "continuous", sampler = "laplace", cellsize = 1.5,
                control.mcmc = ctrl)
  expect_true(nrow(pc$pred.loc) > 0)
  expect_true(all(pc$RR_mean > 0))

  ex <- exceedance(pd, thresholds = c(1, 1.5))
  expect_equal(dim(ex), c(sim$N, 2))
  expect_true(all(ex >= 0 & ex <= 1))
})
