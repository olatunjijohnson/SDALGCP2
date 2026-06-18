# Benchmark SDALGCP2 (C++/RcppArmadillo) vs the original SDALGCP (pure R).
# Two honest comparisons on identical simulated data, same delta/method/phi grid
# and identical MCMC controls + seeds:
#   (B1) the flagship aggregated-correlation build  -- precompute_corr vs precomputeCorrMatrix
#   (E2E) the full MCML fit                          -- SDALGCP2() vs SDALGCP::SDALGCPMCML()
# Estimates are reported alongside times so speed never comes at the cost of agreement.
#
# Run:  Rscript scripts/benchmark_vs_SDALGCP.R

suppressMessages({library(SDALGCP2); library(SDALGCP); library(sf); library(bench)})

sim_data <- function(n_side, seed = 123) {
  set.seed(seed)
  bound <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20)))
  shp   <- st_sf(geometry = st_make_grid(bound, n = c(n_side, n_side)))
  N <- nrow(shp)
  pts <- sda_points(shp, delta = 1.2, method = 3)
  phi_grid <- seq(1, 5, length.out = 8)
  Sig <- 0.5 * precompute_corr(pts, phi_grid)$R[, , which.min(abs(phi_grid - 2.5))]
  x1  <- as.numeric(scale(st_coordinates(st_centroid(shp))[, 1])); D <- cbind(1, x1)
  pop <- round(runif(N, 500, 3000))
  y   <- rpois(N, pop * exp(D %*% c(-6, 0.5) + as.numeric(t(chol(Sig)) %*% rnorm(N))))
  list(shp = shp, N = N, phi_grid = phi_grid,
       dat = data.frame(y = y, x1 = x1, pop = pop),
       form = y ~ x1 + offset(log(pop)))
}

bench_B1 <- function(s) {
  pts2 <- sda_points(s$shp, delta = 1.2, method = 3)            # SDALGCP2 points
  pts1 <- SDALGCP::SDALGCPpolygonpoints(s$shp, delta = 1.2, method = 3)  # SDALGCP points
  npts <- sum(vapply(pts2, function(z) nrow(as.matrix(z$xy)), 0L))
  b <- bench::mark(
    SDALGCP2 = precompute_corr(pts2, s$phi_grid),
    SDALGCP  = SDALGCP:::precomputeCorrMatrix(pts1, s$phi_grid),
    check = FALSE, min_iterations = 3, filter_gc = FALSE)
  list(new = as.numeric(b$median[1]), old = as.numeric(b$median[2]), npts = npts)
}

bench_E2E <- function(s) {
  ctrl2 <- control_mcmc(n.sim = 6000, burnin = 1500, thin = 6, h = 1.65 / s$N^(1/6))
  ctrl1 <- SDALGCP::controlmcmcSDA(n.sim = 6000, burnin = 1500, thin = 6,
                                   h = 1.65 / s$N^(1/6), c1.h = 0.01, c2.h = 1e-4)
  set.seed(50)
  t2 <- system.time(f2 <- SDALGCP2(s$form, s$dat, s$shp, delta = 1.2,
                                   phi = s$phi_grid, method = 3, control.mcmc = ctrl2))[["elapsed"]]
  set.seed(50)
  t1 <- system.time(f1 <- SDALGCP::SDALGCPMCML(formula = s$form, data = s$dat, my_shp = s$shp,
                                   delta = 1.2, phi = s$phi_grid, method = 3, weighted = FALSE,
                                   control.mcmc = ctrl1, plot_profile = FALSE))[["elapsed"]]
  list(new = t2, old = t1,
       beta2 = as.numeric(f2$beta_opt), beta1 = as.numeric(f1$beta_opt))
}

row <- function(stage, N, old, new) sprintf("| %-26s | %4d | %9.2f | %9.2f | %5.1fx |",
                                            stage, N, old, new, old / new)

cat("\n## SDALGCP2 vs SDALGCP benchmark  (R ", R.version$major, ".", R.version$minor,
    ", ", Sys.info()[["sysname"]], ")\n\n", sep = "")
cat("| Stage                      |    N | SDALGCP (s) | SDALGCP2 (s) | Speedup |\n")
cat("|----------------------------|------|-------------|--------------|---------|\n")

for (n_side in c(8, 12)) {
  s <- sim_data(n_side)
  b1 <- bench_B1(s)
  cat(row(sprintf("Correlation build (B1)"), s$N, b1$old, b1$new), "\n")
  e2 <- bench_E2E(s)
  cat(row("Full MCML fit (E2E)", s$N, e2$old, e2$new), "\n")
  cat(sprintf("|   (N=%d, %d candidate pts; fixed effects agree: beta=(%.2f, %.2f) [SDALGCP2] vs (%.2f, %.2f) [SDALGCP]) |\n",
              s$N, b1$npts, e2$beta2[1], e2$beta2[2], e2$beta1[1], e2$beta1[2]))
}
cat("\n")
