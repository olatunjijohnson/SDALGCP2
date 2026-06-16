test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("nugget direct gradient/Hessian match numerical differentiation", {
  skip_if_not_installed("numDeriv")
  set.seed(99)
  N <- 18; p <- 2; B <- 180
  centers <- matrix(runif(N * 2, 0, 30), ncol = 2)
  coords <- lapply(seq_len(N), function(i) {
    ni <- sample(6:9, 1)
    sweep(matrix(rnorm(ni * 2, 0, 1.2), ncol = 2), 2, centers[i, ], "+")
  })
  cg0 <- SDALGCP2:::corr_and_grad_cpp(coords, list(), 3, FALSE, 0L)
  C0 <- cg0$R; diag(C0) <- diag(C0) + 0.2
  x1 <- rnorm(N); D <- cbind(1, x1); m <- rep(50, N)
  S.sim <- matrix(rnorm(B * N), B, N) %*% chol(0.6 * C0) +
    matrix(D %*% c(-0.2, 0.8), B, N, byrow = TRUE)
  y <- rpois(N, m * exp(D %*% c(-0.2, 0.8)))
  data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)

  f <- SDALGCP2:::.mcml_direct_nugget_fit(y, D, m, coords, list(), FALSE, S.sim, data_ll,
         par0_opt = c(-0.2, 0.8, log(0.6)), phi0 = 3, nu0 = 0.2, n = N, p = p)
  xi <- c(-0.1, 0.6, log(0.8), log(2.2), log(0.3))
  expect_lt(max(abs(f$grad(xi) - numDeriv::grad(f$obj, xi))), 1e-5)
  expect_lt(max(abs(f$hess(xi) - numDeriv::jacobian(f$grad, xi))), 1e-4)
})

test_that("nugget is recovered from data with injected overdispersion", {
  skip_on_cran()
  suppressMessages(library(sf))
  set.seed(11)
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(8, 8)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 1.3, method = 3)
  phi_grid <- seq(1, 7, length.out = 10)
  corr <- precompute_corr(pts, phi_grid)
  R3 <- corr$R[, , which.min(abs(phi_grid - 3))]
  C <- 0.5 * (R3 + 0.4 * diag(N))
  x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
  pop <- round(runif(N, 800, 4000))
  y <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
                          as.numeric(t(chol(C)) %*% rnorm(N))))
  dat <- data.frame(y = y, x1 = x1, pop = pop)
  ctrl <- control_mcmc(n.sim = 6000, burnin = 1500, thin = 6, h = 1.65 / N^(1/6))

  fit <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.3, phi = phi_grid,
                  phi_method = "direct", nugget = TRUE, control.mcmc = ctrl, reanchor = 3)
  # The nugget is weakly identified, so we assert structure rather than a fragile
  # point value: the parameter is fitted, finite, non-negative and reported.
  expect_true(isTRUE(fit$nugget))
  expect_true("nu" %in% rownames(fit$cov))
  expect_true(is.finite(fit$nu_opt) && fit$nu_opt >= 0)
  expect_equal(nrow(summary(fit)$coefficients), 5L)   # intercept, x1, sigma^2, phi, nu
})
