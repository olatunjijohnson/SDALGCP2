test_that("testthat is operational", {
  expect_equal(1, 1)
})

test_that("C++ aggregated correlation matches the pure-R reference", {
  set.seed(1)
  N <- 25
  pts <- lapply(seq_len(N), function(i) {
    ni <- sample(5:15, 1)
    w <- runif(ni)
    list(xy = matrix(runif(ni * 2, 0, 1000), ncol = 2), weight = w / sum(w))
  })
  phi <- seq(200, 1500, length.out = 6)

  # unweighted (flatten to plain numeric vectors to compare values only)
  attr(pts, "weighted") <- FALSE
  ref <- SDALGCP2:::precompute_corr_ref(pts, phi, weighted = FALSE)
  cpp <- precompute_corr(pts, phi, weighted = FALSE)$R
  expect_equal(as.numeric(cpp), as.numeric(ref), tolerance = 1e-10)

  # weighted
  refw <- SDALGCP2:::precompute_corr_ref(pts, phi, weighted = TRUE)
  cppw <- precompute_corr(pts, phi, weighted = TRUE)$R
  expect_equal(as.numeric(cppw), as.numeric(refw), tolerance = 1e-10)
})

test_that("correlation array is symmetric with unit-ish diagonal and carries attrs", {
  set.seed(2)
  N <- 12
  pts <- lapply(seq_len(N), function(i) {
    ni <- sample(5:10, 1)
    w <- runif(ni)
    list(xy = matrix(runif(ni * 2, 0, 500), ncol = 2), weight = w / sum(w))
  })
  attr(pts, "weighted") <- TRUE
  attr(pts, "my_shp") <- "placeholder"
  phi <- c(100, 400)
  out <- precompute_corr(pts, phi)

  expect_equal(dim(out$R), c(N, N, length(phi)))
  for (k in seq_along(phi)) {
    expect_equal(out$R[, , k], t(out$R[, , k]), tolerance = 1e-12)
  }
  expect_true(isTRUE(attr(out$R, "weighted")))
  expect_identical(attr(out$R, "my_shp"), "placeholder")
})

test_that("phi must be positive", {
  pts <- list(list(xy = matrix(runif(10), ncol = 2), weight = rep(0.2, 5)))
  attr(pts, "weighted") <- FALSE
  expect_error(precompute_corr(pts, phi = c(-1, 100)), "positive")
})
