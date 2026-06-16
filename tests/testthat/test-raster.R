test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("tilting reduces to areal averaging at beta = 0", {
  set.seed(1)
  Zlist <- list(cbind(`(Intercept)` = 1, z = c(0.2, 0.8, 1.5)),
                cbind(`(Intercept)` = 1, z = c(-1, 0, 1)))
  wlist <- list(c(0.2, 0.5, 0.3), c(1, 1, 1) / 3)
  tl <- SDALGCP2:::.tilt(Zlist, wlist, beta = c(0, 0))
  # b_i = log(sum w) = 0; effective design = areal (w-weighted) mean
  expect_equal(tl$b, c(0, 0), tolerance = 1e-12)
  expect_equal(tl$Dc[, "z"],
               c(sum(wlist[[1]] * Zlist[[1]][, "z"]),
                 sum(wlist[[2]] * Zlist[[2]][, "z"])), tolerance = 1e-12)
})

test_that("SDALGCP2_raster recovers beta where naive areal averaging is biased", {
  skip_on_cran()
  suppressMessages({library(sf); library(terra)})
  set.seed(3)
  r <- rast(xmin = 0, xmax = 18, ymin = 0, ymax = 18, resolution = 0.1)
  xy <- xyFromCell(r, 1:ncell(r))
  src <- matrix(runif(2 * 10, 1, 17), ncol = 2)
  values(r) <- rowSums(sapply(seq_len(nrow(src)), function(s)
    3.0 * exp(-((xy[, 1] - src[s, 1])^2 + (xy[, 2] - src[s, 2])^2) / 0.5)))
  names(r) <- "z"
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 18, ymax = 18))), n = c(6, 6)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 0.6, method = 3)
  Z <- lapply(pts, function(p) cbind(1, terra::extract(r, as.matrix(p$xy))[, "z"]))
  w <- lapply(pts, function(p) p$weight)
  b_true <- sapply(seq_len(N), function(i) log(sum(w[[i]] * exp(as.numeric(Z[[i]] %*% c(-6, 1))))))
  pop <- round(runif(N, 3000, 7000))
  y <- rpois(N, pop * exp(b_true))                 # no spatial noise: isolate covariate
  zbar <- sapply(seq_len(N), function(i) sum(w[[i]] * Z[[i]][, 2]))
  dat <- data.frame(y = y, z = zbar, pop = pop)

  naive <- stats::glm(y ~ z + offset(log(pop)), poisson, dat)
  ctrl <- control_mcmc(n.sim = 3000, burnin = 800, thin = 5, h = 1.65 / N^(1/6))
  fit <- SDALGCP2_raster(y ~ z + offset(log(pop)), dat, shp, delta = 0.6, rasters = r,
                         phi = seq(1.5, 5, length.out = 5), control.mcmc = ctrl, max_iter = 8)

  expect_true(abs(coef(naive)["z"] - 1) > 0.3)              # naive badly biased
  expect_true(abs(fit$beta_opt["z"] - 1) < abs(coef(naive)["z"] - 1) / 2)  # >=2x better
  expect_true("z" %in% names(fit$beta_opt))
  expect_true(isTRUE(fit$raster))

  # fully covariate-tilted variant runs and stays close to the offset-only fit
  fit_t <- SDALGCP2_raster(y ~ z + offset(log(pop)), dat, shp, delta = 0.6, rasters = r,
                           phi = seq(1.5, 5, length.out = 5), control.mcmc = ctrl,
                           tilt_spatial = TRUE, max_iter = 4)
  expect_true(is.finite(fit_t$beta_opt["z"]))
  expect_lt(abs(fit_t$beta_opt["z"] - fit$beta_opt["z"]), 0.5)
})
