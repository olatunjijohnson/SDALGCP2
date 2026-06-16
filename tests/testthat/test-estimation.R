test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("Laplace Newton mode is a stationary point of the conditional density", {
  set.seed(1)
  n <- 30
  loc <- matrix(runif(n * 2, 0, 10), ncol = 2)
  Sigma <- 0.8 * exp(-as.matrix(dist(loc)) / 2)
  mu <- rep(0.1, n); m <- rep(5, n)
  Strue <- as.numeric(t(chol(Sigma)) %*% rnorm(n)) + mu
  y <- rpois(n, m * exp(Strue))

  lap <- SDALGCP2:::laplace_mode_poisson_cpp(y, m, mu, Sigma)
  # gradient of f(S) = sum(y S - m e^S) - 0.5 (S-mu)' Sigma^{-1}(S-mu) at the mode
  g <- (y - m * exp(lap$mode)) - solve(Sigma, lap$mode - mu)
  expect_lt(max(abs(g)), 1e-6)
  # Sigma.tilde = (Sigma^{-1} + diag(m e^mode))^{-1}
  St <- solve(solve(Sigma) + diag(as.numeric(m * exp(lap$mode))))
  expect_equal(as.numeric(lap$Sigma_tilde), as.numeric(St), tolerance = 1e-8)
})

test_that("vectorised MC log-likelihood matches the naive textbook form", {
  set.seed(3)
  n <- 20; p <- 2; nsim <- 150
  D <- cbind(1, runif(n))
  R <- exp(-as.matrix(dist(matrix(runif(n * 2, 0, 8), ncol = 2))) / 2)
  m <- rep(4, n)
  S.sim <- matrix(rnorm(nsim * n), nsim, n) %*% chol(0.7 * R)
  y <- rpois(n, m * exp(as.numeric(D %*% c(0.2, -0.4))))
  data_ll <- as.numeric(S.sim %*% y) - as.numeric(exp(S.sim) %*% m)
  ch <- chol(R); Rinv <- chol2inv(ch); ldetR <- 2 * sum(log(diag(ch)))
  par <- c(0.1, -0.3, log(0.9))
  Den <- SDALGCP2:::.mc_num_loglik(par, D, Rinv, ldetR, S.sim, data_ll, n, p)$num

  obj_v <- function(pp) {
    cp <- SDALGCP2:::.mc_num_loglik(pp, D, Rinv, ldetR, S.sim, data_ll, n, p)
    log(mean(exp(cp$num - Den)))
  }
  naive <- function(pp) {
    beta <- pp[1:p]; mu <- D %*% beta; s2 <- exp(pp[p + 1])
    num <- vapply(seq_len(nsim), function(i) {
      S <- S.sim[i, ]; dS <- S - mu
      sum(y * S - m * exp(S)) -
        0.5 * (n * log(s2) + ldetR + as.numeric(t(dS) %*% Rinv %*% dS) / s2)
    }, numeric(1))
    log(mean(exp(num - Den)))
  }
  pp <- c(-0.1, 0.5, log(1.2))
  expect_equal(obj_v(pp), naive(pp), tolerance = 1e-10)
})

test_that("mcml_fit recovers parameters on simulated areal data", {
  set.seed(2024)
  N <- 50
  centers <- matrix(runif(N * 2, 0, 30), ncol = 2)
  points <- lapply(seq_len(N), function(i) {
    ni <- sample(6:10, 1)
    list(xy = sweep(matrix(rnorm(ni * 2, 0, 1.2), ncol = 2), 2, centers[i, ], "+"),
         weight = rep(1 / ni, ni))
  })
  attr(points, "weighted") <- FALSE
  phi_grid <- seq(1.5, 6, length.out = 8)
  corr <- precompute_corr(points, phi_grid)

  k <- which.min(abs(phi_grid - 3))
  Sig <- 0.6 * corr$R[, , k]
  x1 <- rnorm(N); D <- cbind(1, x1)
  m <- round(runif(N, 80, 200))
  y <- rpois(N, m * exp(D %*% c(-0.2, 0.8) + as.numeric(t(chol(Sig)) %*% rnorm(N))))
  dat <- data.frame(y = y, x1 = x1, pop = m)
  ctrl <- control_mcmc(n.sim = 6000, burnin = 1500, thin = 6, h = 1.65 / N^(1 / 6))

  fit <- mcml_fit(y ~ x1 + offset(log(pop)), dat, corr, control.mcmc = ctrl)
  expect_s3_class(fit, "SDALGCP2")
  # intercept and slope should be in the right ballpark
  expect_lt(abs(fit$beta_opt[2] - 0.8), 0.25)
  expect_true(fit$sigma2_opt > 0)
  sm <- summary(fit)
  expect_equal(rownames(sm$coefficients), c("(Intercept)", "x1", "sigma^2"))
})
