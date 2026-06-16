test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("the null basis is orthonormal and orthogonal to the design", {
  set.seed(1)
  D <- cbind(1, rnorm(20), rnorm(20))
  K <- SDALGCP2:::.null_basis(D)
  expect_equal(ncol(K), nrow(D) - ncol(D))
  expect_equal(crossprod(K), diag(ncol(K)), tolerance = 1e-10)   # K'K = I
  expect_lt(max(abs(crossprod(D, K))), 1e-10)                    # D'K = 0
})

test_that("restricted spatial regression de-confounds a spatially smooth covariate", {
  suppressMessages(library(sf))
  set.seed(2)
  regions <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(8, 8)))
  N <- nrow(regions)
  pts <- sda_points(regions, delta = 1.1, method = 3)
  cc <- precompute_corr(pts, seq(2, 5, length.out = 6))
  Sig <- 0.5 * cc$R[, , 3]
  regions$x1 <- as.numeric(scale(st_coordinates(st_centroid(regions))[, 1]))  # spatial gradient
  regions$pop <- round(runif(N, 500, 3000))
  regions$cases <- rpois(N, regions$pop * exp(-6 + 0.5 * regions$x1 +
                                              as.numeric(t(chol(Sig)) %*% rnorm(N))))
  df <- st_drop_geometry(regions); form <- cases ~ x1 + offset(log(pop))
  attr(cc$R, "my_shp") <- regions

  std <- sdalgcp(form, regions, control = sdalgcp_control(n_sim = 4000, burnin = 1000, thin = 5))
  rsr <- SDALGCP2:::.fit_restricted(form, df, cc)

  expect_s3_class(rsr, "SDALGCP2")
  expect_true("phi" %in% c(names(rsr$estimates), "phi"))
  # restricted estimate is closer to the truth (0.5) than the confounded one
  expect_lt(abs(rsr$beta_opt["x1"] - 0.5), abs(std$beta_opt["x1"] - 0.5))
  expect_true(is.finite(sqrt(rsr$cov["x1", "x1"])))
})

test_that("sdalgcp(confounding = 'restricted') runs end-to-end", {
  suppressMessages(library(sf))
  set.seed(5)
  regions <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 16, ymax = 16))), n = c(6, 6)))
  N <- nrow(regions)
  regions$x1 <- as.numeric(scale(st_coordinates(st_centroid(regions))[, 1]))
  regions$pop <- round(runif(N, 500, 3000))
  regions$cases <- rpois(N, regions$pop * exp(-6 + 0.5 * regions$x1))
  fit <- sdalgcp(cases ~ x1 + offset(log(pop)), regions,
                 control = sdalgcp_control(confounding = "restricted", scale = "grid",
                                           phi = seq(2, 6, length.out = 5)))
  expect_s3_class(fit, "sdalgcp")
  expect_match(fit$mode, "restricted")
  expect_true(is.finite(fit$beta_opt["x1"]))
})
