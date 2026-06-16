suppressMessages({library(SDALGCP2); library(sf)})
set.seed(123)
## simulated lattice of polygons
bound <- st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20)))
grid  <- st_make_grid(bound, n=c(8,8))
shp   <- st_sf(geometry=grid)
N <- nrow(shp); cat("N regions:", N, "\n")

## build candidate points + correlation to simulate a true field
pts  <- sda_points(shp, delta=1.2, method=3)              # regular grid points
phi_grid <- seq(1, 5, length.out=8)
corr <- precompute_corr(pts, phi_grid)
k <- which.min(abs(phi_grid-2.5))
Sig <- 0.5*corr$R[,,k]
x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[,1]))   # a covariate
D  <- cbind(1, x1)
Strue <- as.numeric(t(chol(Sig)) %*% rnorm(N))
pop <- round(runif(N, 500, 3000))
y <- rpois(N, pop*exp(D %*% c(-6, 0.5) + Strue))
dat <- data.frame(y=y, x1=x1, pop=pop)

## full one-call fit
ctrl <- control_mcmc(n.sim=6000, burnin=1500, thin=6, h=1.65/N^(1/6))
t <- system.time(fit <- SDALGCP2(y~x1+offset(log(pop)), dat, shp, delta=1.2,
                  phi=phi_grid, method=3, control.mcmc=ctrl))[["elapsed"]]
cat(sprintf("SDALGCP2() full pipeline: %.2fs\n", t))
print(round(c(fit$beta_opt, sigma2=fit$sigma2_opt, phi=fit$phi_opt),3))

## prediction: discrete (mcmc) + continuous (laplace fast path)
pd <- predict(fit, type="discrete", sampler="mcmc", control.mcmc=ctrl)
cat("discrete RR  range:", round(range(pd$pMean_RR),3), "\n")
pc <- predict(fit, type="continuous", sampler="laplace", cellsize=1.0, control.mcmc=ctrl)
cat("continuous grid pts:", nrow(pc$pred.loc), " RRmean range:", round(range(pc$RRmean),3), "\n")
ex <- exceedance(pd, thresholds=c(1,1.5))
cat("exceedance P(ARR>1) range:", round(range(ex[,1]),3), "\n")
cat("OK\n")
