suppressMessages({library(SDALGCP2); library(sf)})
set.seed(2025)
shp <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(9,9)))
N <- nrow(shp)
pts <- sda_points(shp, delta=1.1, method=3)
phi_grid <- seq(1,7,length.out=15)
corr <- precompute_corr(pts, phi_grid)
k <- which.min(abs(phi_grid-3)); Sig <- 0.5*corr$R[,,k]
cxy <- st_coordinates(st_centroid(shp)); x1 <- as.numeric(scale(cxy[,1]))
pop <- round(runif(N,500,3000))
y <- rpois(N, pop*exp(cbind(1,x1)%*%c(-6,0.5) + as.numeric(t(chol(Sig))%*%rnorm(N))))
dat <- data.frame(y=y,x1=x1,pop=pop)
FORM <- y~x1+offset(log(pop))
ctrl <- control_mcmc(n.sim=8000, burnin=2000, thin=6, h=1.65/N^(1/6))

set.seed(7); tg <- system.time(fg <- SDALGCP2(FORM,dat,shp,delta=1.1,phi=phi_grid,method=3,control.mcmc=ctrl,phi_method="grid"))[["elapsed"]]
set.seed(7); td <- system.time(fd <- SDALGCP2(FORM,dat,shp,delta=1.1,phi=phi_grid,method=3,control.mcmc=ctrl,phi_method="direct"))[["elapsed"]]

cat(sprintf("True:   beta=(-6, 0.5) sigma2=0.5 phi=3\n"))
cat(sprintf("GRID:   beta=(%.3f, %.3f) sigma2=%.3f phi=%.3f   [%.2fs]\n", fg$beta_opt[1],fg$beta_opt[2],fg$sigma2_opt,fg$phi_opt,tg))
cat(sprintf("DIRECT: beta=(%.3f, %.3f) sigma2=%.3f phi=%.3f   [%.2fs]\n", fd$beta_opt[1],fd$beta_opt[2],fd$sigma2_opt,fd$phi_opt,td))
cat(sprintf("phi SE (direct, from Hessian): %.3f\n", sqrt(fd$cov["phi","phi"])))
cat("\n--- summary(direct) includes phi row ---\n"); print(round(summary(fd)$coefficients,4))
