## Tutorial 6 (misaligned covariates) — point + areal, figures + output.
suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})
th <- theme_minimal(base_size=11) + theme(axis.text=element_blank(), panel.grid=element_blank())
zf <- function(xy) 1.5*sin(xy[,1]/4) + 1.2*cos(xy[,2]/5) + 0.5*xy[,1]/20
## outcome regions + true counts from the point-level covariate (true beta=0.8)
set.seed(123)
reg <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(8,8)))
N <- nrow(reg); pts <- sda_points(reg, delta=1.0, method=3); w <- lapply(pts, function(p)p$weight)
b_true <- sapply(1:N, function(i){Z<-cbind(1,zf(as.matrix(pts[[i]]$xy))); log(sum(w[[i]]*exp(as.numeric(Z%*%c(-6,0.8)))))})
reg$pop <- round(runif(N,800,5000)); reg$cases <- rpois(N, reg$pop*exp(b_true))
## POINT covariate: 40 monitors
mon <- st_as_sf(data.frame(x=runif(40,0,20), y=runif(40,0,20)), coords=c("x","y")); mon$z <- zf(st_coordinates(mon))
## kriging-surface figure
gp <- SDALGCP2:::.gp_point(st_coordinates(mon), mon$z)
grid <- expand.grid(x=seq(0.5,19.5,0.5), y=seq(0.5,19.5,0.5)); grid$z <- SDALGCP2:::.krige_point(gp, as.matrix(grid))$mean
mc <- st_coordinates(mon)
g <- ggplot() + geom_raster(data=grid, aes(x,y,fill=z)) +
  geom_point(data=data.frame(x=mc[,1],y=mc[,2],z=mon$z), aes(x,y,fill=z), shape=21, size=3, color="white", stroke=0.6) +
  geom_sf(data=st_union(reg), fill=NA, color="black") + scale_fill_viridis_c(name="covariate") +
  coord_sf(expand=FALSE) + th + labs(title="Covariate kriged from 40 monitors to the outcome's support")
ggsave("vignettes/t6_krige.png", g, width=5.6, height=4.4, dpi=120)
## fits
fit_pt <- sdalgcp(cases ~ z + offset(log(pop)), reg, covariates=list(z=mon))
## AREAL covariate: averages over a different 4x4 partition
covpoly <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(4,4)))
cpts <- sda_points(covpoly, delta=1.0, method=3); covpoly$z <- sapply(cpts, function(p) mean(zf(as.matrix(p$xy))))
fit_ar <- sdalgcp(cases ~ z + offset(log(pop)), reg, covariates=list(z=covpoly))
sink("vignettes/t6_out.txt")
cat(sprintf("True z effect: 0.80\n"))
cat(sprintf("Point support  (40 monitors):           %.3f\n", fit_pt$beta_opt["z"]))
cat(sprintf("Areal support  (25 different polygons):  %.3f\n", fit_ar$beta_opt["z"]))
sink()
cat("T6 done\n")
