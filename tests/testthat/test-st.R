test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("Kronecker-free identities match the brute-force Kronecker product", {
  set.seed(1)
  N <- 6; T <- 4
  Rs <- exp(-as.matrix(dist(matrix(runif(N * 2), ncol = 2))) / 0.3)
  Rt <- SDALGCP2:::.temporal_corr(1:T, 2, 0.5)
  M <- matrix(rnorm(N * T), N, T)
  x <- as.vector(M)
  Big <- kronecker(Rt, Rs)
  q_kf <- sum((solve(Rs) %*% M %*% solve(Rt)) * M)
  ld_kf <- N * as.numeric(determinant(Rt, log = TRUE)$modulus) +
           T * as.numeric(determinant(Rs, log = TRUE)$modulus)
  expect_equal(q_kf, as.numeric(t(x) %*% solve(Big) %*% x), tolerance = 1e-9)
  expect_equal(ld_kf, as.numeric(determinant(Big, log = TRUE)$modulus), tolerance = 1e-9)
})

test_that("temporal_corr derivative matches finite differences", {
  for (kt in c(0.5, 1.5, 2.5)) {
    tc <- SDALGCP2:::.temporal_corr(1:5, 2, kt, deriv = TRUE)
    h <- 1e-5
    fd <- (SDALGCP2:::.temporal_corr(1:5, 2 + h, kt) -
           SDALGCP2:::.temporal_corr(1:5, 2 - h, kt)) / (2 * h)
    expect_lt(max(abs(tc$dR - fd)), 1e-6)
  }
})

test_that("ST analytic gradient matches numerical differentiation", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  N <- 8; T <- 4; p <- 2; B <- 120; times <- 1:T; kt <- 0.5
  coords <- lapply(seq_len(N), function(i) {
    ni <- sample(5:8, 1); matrix(runif(ni * 2, 0, 50), ncol = 2)
  })
  pp <- lapply(coords, function(z) list(xy = z, weight = rep(1 / nrow(z), nrow(z))))
  attr(pp, "weighted") <- FALSE
  Rs <- precompute_corr(pp, c(10, 20))$R[, , 1]
  Rsi <- solve(Rs); ldetRs <- as.numeric(determinant(Rs, log = TRUE)$modulus)
  Rt <- SDALGCP2:::.temporal_corr(times, 2, kt)
  x1 <- rnorm(N * T); D <- cbind(1, x1)
  L <- t(chol(0.6 * kronecker(Rt, Rs)))
  S.sim <- t(apply(matrix(rnorm(B * N * T), B, N * T), 1,
                   function(z) as.numeric(D %*% c(-0.5, 0.3) + L %*% z)))
  y <- rpois(N * T, 5 * exp(as.numeric(D %*% c(-0.5, 0.3)))); m <- rep(5, N * T)
  data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)
  Den <- SDALGCP2:::.st_num_loglik(c(-0.5, 0.3, log(0.6), log(2)), D, Rsi, ldetRs,
                                   S.sim, data_ll, N, T, p, times, kt)$num
  obj <- function(th) -log(mean(exp(
    SDALGCP2:::.st_num_loglik(th, D, Rsi, ldetRs, S.sim, data_ll, N, T, p, times, kt)$num - Den)))
  negGrad <- function(theta) {
    cp <- SDALGCP2:::.st_num_loglik(theta, D, Rsi, ldetRs, S.sim, data_ll, N, T, p,
                                    times, kt, want_grad = TRUE)
    w <- exp(cp$num - Den); w <- w / sum(w); s2 <- cp$s2; nu <- cp$nu
    gb <- t(sapply(seq_along(w), function(b) as.numeric(crossprod(D, as.vector(cp$G[[b]]))) / s2))
    gs <- -N * T / 2 + cp$q / (2 * s2); gnu <- nu * (-0.5 * N * cp$trRtiD + 0.5 * cp$qnu / s2)
    -colSums(w * cbind(gb, gs, gnu))
  }
  th <- c(-0.3, 0.2, log(0.8), log(2.5))
  expect_lt(max(abs(negGrad(th) - numDeriv::grad(obj, th))), 1e-5)
})

test_that("SDALGCP2_ST fits end-to-end and recovers the slope", {
  skip_on_cran()
  suppressMessages(library(sf))
  set.seed(7)
  shp <- st_sf(geometry = st_make_grid(
    st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 15, ymax = 15))), n = c(5, 5)))
  N <- nrow(shp); T <- 4; times <- 1:T
  pts <- sda_points(shp, delta = 1.5, method = 3)
  phi_grid <- seq(1.5, 5, length.out = 6); corr <- precompute_corr(pts, phi_grid)
  Rs <- corr$R[, , 3]; Rt <- SDALGCP2:::.temporal_corr(times, 2, 0.5)
  L <- t(chol(0.4 * kronecker(Rt, Rs)))
  x1 <- as.numeric(scale(rep(st_coordinates(st_centroid(shp))[, 1], T)))
  pop <- round(runif(N * T, 800, 4000))
  y <- rpois(N * T, pop * exp(cbind(1, x1) %*% c(-6, 0.5) + as.numeric(L %*% rnorm(N * T))))
  dat <- data.frame(y = y, x1 = x1, pop = pop)
  ctrl <- control_mcmc(n.sim = 4000, burnin = 1000, thin = 5, h = 1.65 / (N * T)^(1/6))
  fit <- SDALGCP2_ST(y ~ x1 + offset(log(pop)), dat, shp, times = times, delta = 1.5,
                     phi = phi_grid, control.mcmc = ctrl)
  expect_s3_class(fit, "SDALGCP2_ST")
  expect_lt(abs(fit$beta_opt[2] - 0.5), 0.3)
  expect_true(is.finite(fit$nu_opt) && fit$nu_opt > 0)
  expect_true("nu" %in% rownames(fit$cov))

  # re-anchoring runs and discrete prediction returns N x T risk surfaces
  fit2 <- SDALGCP2_ST(y ~ x1 + offset(log(pop)), dat, shp, times = times, delta = 1.5,
                      phi = phi_grid, control.mcmc = ctrl, reanchor = 2)
  expect_true(fit2$n_reanchor >= 1)
  pred <- predict(fit2)
  # prediction is a long sf (region x time), mirroring the spatial predictor
  expect_s3_class(pred, "SDALGCP2_ST_pred")
  expect_s3_class(pred, "sf")
  expect_equal(nrow(pred), N * T)
  expect_true(all(c("region", "time", "relative_risk", "adjusted_rr") %in% names(pred)))
  expect_true(all(pred$relative_risk > 0))
  expect_true(all(pred$adjusted_rr > 0))
  expect_equal(attr(pred, "N"), N)
  expect_equal(attr(pred, "T"), T)
  # a time slice maps without error
  expect_s3_class(plot(pred, time = times[1]), "ggplot")
})
