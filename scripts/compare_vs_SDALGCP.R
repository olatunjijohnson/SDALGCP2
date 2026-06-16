suppressMessages({library(SDALGCP2); library(SDALGCP); library(sf)})
set.seed(123)
bound <- st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20)))
shp   <- st_sf(geometry=st_make_grid(bound, n=c(8,8)))
N <- nrow(shp)
pts  <- sda_points(shp, delta=1.2, method=3)
phi_grid <- seq(1, 5, length.out=8)
corr <- precompute_corr(pts, phi_grid)
k <- which.min(abs(phi_grid-2.5)); Sig <- 0.5*corr$R[,,k]
x1 <- as.numeric(scale(st_coordinates(st_centroid(shp))[,1])); D <- cbind(1,x1)
pop <- round(runif(N,500,3000))
y <- rpois(N, pop*exp(D%*%c(-6,0.5) + as.numeric(t(chol(Sig))%*%rnorm(N))))
dat <- data.frame(y=y, x1=x1, pop=pop)
FORM <- y~x1+offset(log(pop))
ctrl2 <- control_mcmc(n.sim=6000, burnin=1500, thin=6, h=1.65/N^(1/6))
ctrl1 <- SDALGCP::controlmcmcSDA(n.sim=6000, burnin=1500, thin=6, h=1.65/N^(1/6), c1.h=0.01, c2.h=1e-4)

set.seed(50)
t2 <- system.time(f2 <- SDALGCP2(FORM, dat, shp, delta=1.2, phi=phi_grid, method=3,
                    control.mcmc=ctrl2))[["elapsed"]]
set.seed(50)
t1 <- system.time(f1 <- SDALGCP::SDALGCPMCML(formula=FORM, data=dat, my_shp=shp, delta=1.2,
                    phi=phi_grid, method=3, weighted=FALSE, control.mcmc=ctrl1,
                    plot_profile=FALSE))[["elapsed"]]

fmt <- function(b,s,p) paste0("beta=(",paste(round(b,3),collapse=", "),") sigma2=",round(s,3)," phi=",round(p,3))
cat("True:     ", fmt(c(-6,0.5),0.5,2.5), "\n")
cat("SDALGCP2: ", fmt(f2$beta_opt, f2$sigma2_opt, f2$phi_opt), sprintf("  [%.2fs]\n",t2))
cat("SDALGCP:  ", fmt(as.numeric(f1$beta_opt), f1$sigma2_opt, f1$phi_opt), sprintf("  [%.2fs]\n",t1))
cat(sprintf("\nFull-pipeline speedup: %.1fx\n", t1/t2))
