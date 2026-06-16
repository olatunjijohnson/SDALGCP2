test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("sdalgcp_control defaults and auto-delta are sensible", {
  ctl <- sdalgcp_control()
  expect_equal(ctl$scale, "continuous")
  expect_equal(ctl$point_method, "regular")
  expect_equal(ctl$reanchor, 2L)
  suppressMessages(library(sf))
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(8, 8)))
  d <- SDALGCP2:::.auto_delta(shp, 16)
  expect_true(d > 0 && is.finite(d))
})

test_that("the minimal sdalgcp() call fits and predict() returns an sf", {
  suppressMessages(library(sf))
  set.seed(2)
  regions <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 18, ymax = 18))), n = c(6, 6)))
  N <- nrow(regions)
  pts <- sda_points(regions, delta = 1.4, method = 3); cc <- precompute_corr(pts, c(2, 3, 4))
  Sig <- 0.5 * cc$R[, , 2]
  regions$x1 <- as.numeric(scale(st_coordinates(st_centroid(regions))[, 1]))
  regions$pop <- round(runif(N, 500, 3000))
  regions$cases <- rpois(N, regions$pop * exp(-6 + 0.5 * regions$x1 +
                                              as.numeric(t(chol(Sig)) %*% rnorm(N))))

  fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = regions,
                 control = sdalgcp_control(n_sim = 4000, burnin = 1000, thin = 5))
  expect_s3_class(fit, "sdalgcp")
  expect_true(is.finite(fit$phi_opt) && fit$phi_opt > 0)
  expect_output(print(fit), "Spatially discrete LGCP")

  rr <- predict(fit)
  expect_s3_class(rr, "sf")
  expect_true(all(c("relative_risk", "relative_risk_se", "incidence") %in% names(rr)))
  expect_true(all(rr$relative_risk > 0))

  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})

test_that("sdalgcp() rejects non-sf data with a helpful message", {
  expect_error(sdalgcp(y ~ x, data = data.frame(y = 1:3, x = 1:3)), "sf")
})
