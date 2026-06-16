# Spatially continuous (raster) covariates: why intensity-scale aggregation beats
# averaging the predictor over polygons.
#
# We simulate a covariate with sharp within-region "exposure source" peaks (e.g.
# pollution sources). A region containing a peak has a modest AREAL-MEAN covariate
# but a large intensity-scale aggregate sum_k w_ik exp(beta z_ik). Fitting a
# Poisson model on the areal-mean covariate is therefore badly biased, while
# SDALGCP2_raster() -- which aggregates on the intensity scale -- recovers beta.

suppressMessages({library(sf); library(terra); library(SDALGCP2); library(ggplot2)})
set.seed(3)

## ---- a covariate with sharp sub-region peaks ----
r <- rast(xmin = 0, xmax = 20, ymin = 0, ymax = 20, resolution = 0.08)
xy <- xyFromCell(r, 1:ncell(r))
src <- matrix(runif(2 * 14, 1, 19), ncol = 2)
values(r) <- rowSums(sapply(seq_len(nrow(src)), function(s)
  3.2 * exp(-((xy[, 1] - src[s, 1])^2 + (xy[, 2] - src[s, 2])^2) / 0.5)))
names(r) <- "z"

shp <- st_sf(geometry = st_make_grid(
  st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 20, ymax = 20))), n = c(8, 8)))
N <- nrow(shp)

## ---- simulate aggregated counts from the point-level intensity model ----
pts <- sda_points(shp, delta = 0.5, method = 3)
Z   <- lapply(pts, function(p) cbind(1, terra::extract(r, as.matrix(p$xy))[, "z"]))
w   <- lapply(pts, function(p) p$weight)
beta_true <- c(-6, 1.0)
b_true <- sapply(seq_len(N), function(i) log(sum(w[[i]] * exp(as.numeric(Z[[i]] %*% beta_true)))))
pop <- round(runif(N, 2000, 6000))
y   <- rpois(N, pop * exp(b_true))
zbar <- sapply(seq_len(N), function(i) sum(w[[i]] * Z[[i]][, 2]))   # areal-mean covariate
dat <- data.frame(y = y, z = zbar, pop = pop)

## ---- (a) naive areal-average covariate ; (b) intensity-scale ----
naive <- glm(y ~ z + offset(log(pop)), poisson, dat)
ctrl  <- control_mcmc(n.sim = 6000, burnin = 1500, thin = 6, h = 1.65 / N^(1/6))
fit_r <- SDALGCP2_raster(y ~ z + offset(log(pop)), dat, shp, delta = 0.5, rasters = r,
                         phi = seq(1.5, 6, length.out = 8), control.mcmc = ctrl, max_iter = 10)

cat(sprintf("\nTrue z effect:        1.000\n"))
cat(sprintf("NAIVE areal-average:  %.3f   (bias %+.0f%%)\n", coef(naive)["z"], 100 * (coef(naive)["z"] - 1)))
cat(sprintf("SDALGCP2_raster:      %.3f   (bias %+.0f%%)\n", fit_r$beta_opt["z"], 100 * (fit_r$beta_opt["z"] - 1)))

## ---- figure: the covariate surface with region outlines ----
rdf <- as.data.frame(r, xy = TRUE)
g <- ggplot() +
  geom_raster(data = rdf, aes(x, y, fill = z)) +
  geom_sf(data = shp, fill = NA, color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(name = "z(x)") +
  coord_sf(expand = FALSE) + theme_minimal() +
  labs(title = "Raster covariate with sharp sub-region peaks",
       subtitle = "Areal averaging washes out the peaks; intensity-scale aggregation does not")
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("man/figures/raster_covariate.png", g, width = 5.2, height = 4.4, dpi = 110)
