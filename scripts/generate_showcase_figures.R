suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})
set.seed(42)
## 10x10 regions; covariate is NON-spatial (avoids confounding); add a spatial hotspot field
regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(10,10)))
N <- nrow(regions)
cxy <- st_coordinates(st_centroid(regions))
pts <- sda_points(regions, delta=0.9, method=3); cc <- precompute_corr(pts, c(2,3,4,5))
Sig <- 0.5*cc$R[,,3]                                  # true phi ~ 4
S <- as.numeric(t(chol(Sig)) %*% rnorm(N))
regions$x1  <- rnorm(N)                               # non-spatial covariate
regions$pop <- round(runif(N, 800, 5000))
regions$cases <- rpois(N, regions$pop*exp(-6 + 0.6*regions$x1 + S))

fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data = regions,
               control = sdalgcp_control(reanchor=3))
cat(sprintf("beta x1 = %.3f (true 0.6), phi = %.2f (true 4)\n", fit$beta_opt["x1"], fit$phi_opt))

th <- theme_minimal(base_size=12) + theme(axis.text=element_blank(), panel.grid=element_blank())
g1 <- plot(fit, "risk") + th
g2 <- plot(fit, "risk_se") + th
g3 <- plot(fit, "exceedance", threshold=1.5) + th
dir.create("man/figures", showWarnings=FALSE, recursive=TRUE)
ggsave("man/figures/showcase_risk.png", g1, width=5, height=4.2, dpi=120)
ggsave("man/figures/showcase_uncertainty.png", g2, width=5, height=4.2, dpi=120)
ggsave("man/figures/showcase_exceedance.png", g3, width=5, height=4.2, dpi=120)

## continuous surface
pc <- predict(SDALGCP2:::.strip_sdalgcp(fit), type="continuous", sampler="laplace", cellsize=0.6)
g4 <- plot(pc, "relative_risk", bound=regions) + th + labs(title="Continuous relative-risk surface")
ggsave("man/figures/showcase_continuous.png", g4, width=5, height=4.2, dpi=120)
cat("figures written\n")
