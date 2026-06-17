test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("plotting, profile and diagnostics return the right object types", {
  suppressMessages(library(sf))
  set.seed(11)
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 16, ymax = 16))), n = c(6, 6)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 1.5, method = 3)
  phi_grid <- seq(1, 5, length.out = 6)
  corr <- precompute_corr(pts, phi_grid)
  Sig <- 0.5 * corr$R[, , which.min(abs(phi_grid - 2.5))]
  x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
  pop <- round(runif(N, 500, 3000))
  y <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
                          as.numeric(t(chol(Sig)) %*% rnorm(N))))
  dat <- data.frame(y = y, x1 = x1, pop = pop)
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / N^(1/6))
  fit <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.5,
                  phi = phi_grid, method = 3, control.mcmc = ctrl)

  pd <- predict(fit, type = "discrete", control.mcmc = ctrl)
  expect_s3_class(plot(pd, "ARR", midpoint = 1), "ggplot")
  expect_s3_class(plot(pd, "RR_se"), "ggplot")
  expect_s3_class(map_exceedance(pd, threshold = 1.2), "ggplot")
  expect_s3_class(coef_plot(fit), "ggplot")

  pp <- phi_profile(fit, plot = FALSE)
  expect_length(pp$ci, 2)
  expect_true(pp$ci[1] <= fit$phi_opt + 1e-8 && pp$ci[2] >= fit$phi_opt - 1e-8)

  mc <- model_check(fit, pd, nsim = 99, plot = FALSE)
  expect_length(mc$residuals, N)
  expect_true(is.finite(mc$moran$I))

  rp <- report(fit, pd)
  expect_named(rp, c("relative_risk", "uncertainty", "exceedance",
                     "coefficients", "phi_profile"))
})
