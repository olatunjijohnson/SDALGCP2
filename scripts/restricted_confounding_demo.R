suppressMessages({library(SDALGCP2); library(sf)})
set.seed(2)
regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(8,8)))
N <- nrow(regions)
pts <- sda_points(regions, delta=1.1, method=3); cc <- precompute_corr(pts, seq(2,5,length.out=6))
Sig <- 0.5*cc$R[,,3]
## SPATIALLY SMOOTH covariate (a gradient) -> confounding
regions$x1 <- as.numeric(scale(st_coordinates(st_centroid(regions))[,1]))
regions$pop <- round(runif(N,500,3000))
regions$cases <- rpois(N, regions$pop*exp(-6 + 0.5*regions$x1 + as.numeric(t(chol(Sig))%*%rnorm(N))))
df <- st_drop_geometry(regions)
form <- cases ~ x1 + offset(log(pop))
attr(cc$R,"my_shp") <- regions

glm_b <- coef(glm(form, poisson, df))["x1"]
std <- sdalgcp(form, regions)                                  # standard (confounded)
rsr <- SDALGCP2:::.fit_restricted(form, df, cc)               # restricted
cat(sprintf("True x1 = 0.500\n"))
cat(sprintf("Poisson GLM (no spatial):   %.3f\n", glm_b))
cat(sprintf("Standard SDA-LGCP:          %.3f  (confounded, attenuated)\n", std$beta_opt["x1"]))
cat(sprintf("Restricted (RSR):           %.3f  (SE %.3f)\n", rsr$beta_opt["x1"], sqrt(rsr$cov["x1","x1"])))
cat(sprintf("  sigma2=%.3f phi=%.2f\n", rsr$sigma2_opt, rsr$phi_opt))
cat("\n--- via the easy interface ---\n")
rsr2 <- sdalgcp(form, regions, control=sdalgcp_control(confounding="restricted", scale="grid", phi=seq(2,6,length.out=6)))
print(rsr2)
