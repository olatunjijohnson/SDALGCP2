suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})

## ---------- RSR / confounding comparison ----------
set.seed(2)
regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(8,8)))
N <- nrow(regions); pts <- sda_points(regions, delta=1.1, method=3)
cc <- precompute_corr(pts, seq(2,5,length.out=6)); Sig <- 0.5*cc$R[,,3]
regions$x1 <- as.numeric(scale(st_coordinates(st_centroid(regions))[,1]))  # spatial gradient
regions$pop <- round(runif(N,500,3000))
regions$cases <- rpois(N, regions$pop*exp(-6 + 0.5*regions$x1 + as.numeric(t(chol(Sig))%*%rnorm(N))))
df <- st_drop_geometry(regions); attr(cc$R,"my_shp") <- regions
glm_b <- coef(glm(cases~x1+offset(log(pop)), poisson, df))["x1"]
std <- sdalgcp(cases~x1+offset(log(pop)), regions)
rsr <- SDALGCP2:::.fit_restricted(cases~x1+offset(log(pop)), df, cc)
est <- data.frame(method=factor(c("truth","Poisson GLM","standard\nSDA-LGCP","restricted\n(RSR)"),
                  levels=c("truth","Poisson GLM","standard\nSDA-LGCP","restricted\n(RSR)")),
                  beta=c(0.5, glm_b, std$beta_opt["x1"], rsr$beta_opt["x1"]))
g <- ggplot(est, aes(method, beta, fill=method)) + geom_col(width=0.6) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_text(aes(label=sprintf("%.2f",beta)), vjust=-0.4) +
  scale_fill_manual(values=c("grey50","#888888","#B2182B","#2166AC"), guide="none") +
  theme_minimal(base_size=12) + labs(title="Spatial confounding: covariate effect",
    subtitle="Standard fit is attenuated; RSR recovers it", x=NULL, y=expression(hat(beta)[x1]))
ggsave("vignettes/rsr_compare.png", g, width=5.4, height=4, dpi=120)
cat(sprintf("RSR fig: glm=%.2f std=%.2f rsr=%.2f\n", glm_b, std$beta_opt["x1"], rsr$beta_opt["x1"]))

## ---------- Misalignment: kriged covariate surface + monitors ----------
set.seed(123)
zf <- function(xy) 1.5*sin(xy[,1]/4) + 1.2*cos(xy[,2]/5) + 0.5*xy[,1]/20
reg2 <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(8,8)))
mon <- st_as_sf(data.frame(x=runif(40,0,20), y=runif(40,0,20)), coords=c("x","y")); mon$z <- zf(st_coordinates(mon))
gp <- SDALGCP2:::.gp_point(st_coordinates(mon), mon$z)
grid <- expand.grid(x=seq(0.5,19.5,0.5), y=seq(0.5,19.5,0.5))
kr <- SDALGCP2:::.krige_point(gp, as.matrix(grid)); grid$z <- kr$mean
mc <- st_coordinates(mon)
g2 <- ggplot() + geom_raster(data=grid, aes(x,y,fill=z)) +
  geom_point(data=data.frame(x=mc[,1],y=mc[,2],z=mon$z), aes(x,y,fill=z), shape=21, size=3, color="white", stroke=0.6) +
  geom_sf(data=st_union(reg2), fill=NA, color="black") +
  scale_fill_viridis_c(name="covariate") + coord_sf(expand=FALSE) +
  theme_minimal(base_size=12) + theme(axis.text=element_blank(), panel.grid=element_blank()) +
  labs(title="Covariate kriged from monitors to the outcome's support",
       subtitle="Points = 40 monitors; surface = kriged prediction")
ggsave("vignettes/misalign_krige.png", g2, width=5.6, height=4.4, dpi=120)
cat("misalign fig done\n")
