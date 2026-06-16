test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("corr_and_grad_cpp derivatives match finite differences", {
  set.seed(1)
  N <- 15
  coords <- lapply(seq_len(N), function(i) {
    ni <- sample(5:9, 1)
    matrix(runif(ni * 2, 0, 100), ncol = 2)
  })
  phi <- 20; h <- 1e-4
  cg <- SDALGCP2:::corr_and_grad_cpp(coords, list(), phi, FALSE, 0L)
  Rp <- SDALGCP2:::corr_and_grad_cpp(coords, list(), phi + h, FALSE, 0L)$R
  Rm <- SDALGCP2:::corr_and_grad_cpp(coords, list(), phi - h, FALSE, 0L)$R
  expect_lt(max(abs(cg$dR - (Rp - Rm) / (2 * h))), 1e-6)
  expect_lt(max(abs(cg$d2R - (Rp - 2 * cg$R + Rm) / (h * h))), 1e-4)
})

test_that("analytic gradient/Hessian of the direct MCML match numerical", {
  skip_if_not_installed("numDeriv")
  set.seed(99)
  N <- 20; p <- 2; B <- 200
  centers <- matrix(runif(N * 2, 0, 30), ncol = 2)
  pts <- lapply(seq_len(N), function(i) {
    ni <- sample(6:9, 1)
    list(xy = sweep(matrix(rnorm(ni * 2, 0, 1.2), ncol = 2), 2, centers[i, ], "+"),
         weight = rep(1 / ni, ni))
  })
  coords <- lapply(pts, function(z) z$xy)
  cg0 <- SDALGCP2:::corr_and_grad_cpp(coords, list(), 3, FALSE, 0L)
  x1 <- rnorm(N); D <- cbind(1, x1); m <- rep(50, N)
  S.sim <- matrix(rnorm(B * N), B, N) %*% chol(0.6 * cg0$R) +
    matrix(D %*% c(-0.2, 0.8), B, N, byrow = TRUE)
  y <- rpois(N, m * exp(D %*% c(-0.2, 0.8)))
  data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)

  f <- SDALGCP2:::.mcml_direct_fit(y, D, m, coords, list(), FALSE, S.sim, data_ll,
                                   par0_opt = c(-0.2, 0.8, log(0.6)), phi0 = 3,
                                   n = N, p = p)
  xi <- c(-0.1, 0.6, log(0.8), log(2.2))
  expect_lt(max(abs(f$grad(xi) - numDeriv::grad(f$obj, xi))), 1e-5)
  expect_lt(max(abs(f$hess(xi) - numDeriv::jacobian(f$grad, xi))), 1e-4)
})

test_that("direct fit runs end-to-end and gives a phi standard error", {
  suppressMessages(library(sf))
  set.seed(2025)
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 18, ymax = 18))), n = c(7, 7)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 1.3, method = 3)
  phi_grid <- seq(1, 6, length.out = 8)
  corr <- precompute_corr(pts, phi_grid)
  Sig <- 0.5 * corr$R[, , which.min(abs(phi_grid - 3))]
  x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1]))
  pop <- round(runif(N, 500, 3000))
  y <- rpois(N, pop * exp(cbind(1, x1) %*% c(-6, 0.5) +
                          as.numeric(t(chol(Sig)) %*% rnorm(N))))
  dat <- data.frame(y = y, x1 = x1, pop = pop)
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / N^(1/6))

  fit <- SDALGCP2(y ~ x1 + offset(log(pop)), dat, shp, delta = 1.3, phi = phi_grid,
                  method = 3, control.mcmc = ctrl, phi_method = "direct")
  expect_identical(fit$phi_method, "direct")
  expect_true("phi" %in% rownames(fit$cov))
  expect_true(is.finite(sqrt(fit$cov["phi", "phi"])))
  expect_equal(nrow(summary(fit)$coefficients), 4L)   # intercept, x1, sigma^2, phi
})
