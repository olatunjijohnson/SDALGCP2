suppressMessages({library(SDALGCP2); library(sf)})
set.seed(123)
## true covariate field on a fine grid (a GP-like smooth surface)
gx <- expand.grid(x=seq(0,20,0.5), y=seq(0,20,0.5))
zfield <- 1.5*sin(gx$x/4) + 1.2*cos(gx$y/5) + 0.5*gx$x/20
zfun <- function(xy) { # nearest-grid lookup
  i <- apply(xy, 1, function(p) which.min((gx$x-p[1])^2+(gx$y-p[2])^2)); zfield[i] }
## outcome regions
regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(8,8)))
N <- nrow(regions)
pts <- sda_points(regions, delta=1.0, method=3); w <- lapply(pts, function(p)p$weight)
bt <- c(-6, 0.8)   # intercept, z effect
b_true <- sapply(seq_len(N), function(i){ Z <- cbind(1, zfun(as.matrix(pts[[i]]$xy)))
  log(sum(w[[i]]*exp(as.numeric(Z%*%bt)))) })
regions$pop <- round(runif(N,1000,5000)); regions$cases <- rpois(N, regions$pop*exp(b_true))
## covariate observed only at 40 scattered MONITOR points
mon <- st_as_sf(data.frame(x=runif(40,0,20), y=runif(40,0,20)), coords=c("x","y"))
mon$z <- zfun(st_coordinates(mon))
form <- cases ~ z + offset(log(pop))
ctrl <- control_mcmc(n.sim=4000, burnin=1000, thin=5, h=1.65/N^(1/6))
f_naive <- SDALGCP2_misaligned(form, regions, delta=1.0, covariates=list(z=mon), berkson=FALSE, control.mcmc=ctrl, max_iter=8)
f_berk  <- SDALGCP2_misaligned(form, regions, delta=1.0, covariates=list(z=mon), berkson=TRUE,  control.mcmc=ctrl, max_iter=8)
cat(sprintf("True z effect:                 0.800\n"))
cat(sprintf("Naive (kriged-mean plug-in):   %.3f\n", f_naive$beta_opt["z"]))
cat(sprintf("Berkson-corrected:             %.3f\n", f_berk$beta_opt["z"]))
