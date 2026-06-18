## Tutorial 2 (raster covariates) — figures + output.
suppressMessages({library(SDALGCP2); library(sf); library(terra); library(ggplot2)})
th <- theme_minimal(base_size=11) + theme(axis.text=element_blank(), panel.grid=element_blank())
set.seed(3)
## sharp 'exposure source' covariate on a fine raster
r <- rast(xmin=0,xmax=20,ymin=0,ymax=20,resolution=0.08); xy <- xyFromCell(r,1:ncell(r))
src <- matrix(runif(2*14,1,19),ncol=2)
values(r) <- rowSums(sapply(1:nrow(src),function(s)3.2*exp(-((xy[,1]-src[s,1])^2+(xy[,2]-src[s,2])^2)/0.5)))
names(r)<-"z"
sh <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))),n=c(8,8)))
N <- nrow(sh); pp <- sda_points(sh, delta=0.5, method=3); w <- lapply(pp, function(p)p$weight)
Z <- lapply(pp, function(p) cbind(1, terra::extract(r, as.matrix(p$xy))[,"z"]))
bt <- c(-6,1.0)
b_true <- sapply(1:N, function(i) log(sum(w[[i]]*exp(as.numeric(Z[[i]]%*%bt)))))
sh$pop <- round(runif(N,2000,6000)); sh$cases <- rpois(N, sh$pop*exp(b_true))
sh$zbar <- sapply(1:N, function(i) sum(w[[i]]*Z[[i]][,2]))   # areal mean of z

## Fig 1: raster + region borders
rdf <- as.data.frame(r, xy=TRUE)
g1 <- ggplot() + geom_raster(data=rdf, aes(x,y,fill=z)) +
  geom_sf(data=sh, fill=NA, color="white", linewidth=0.25) +
  scale_fill_viridis_c(name="z(x)") + coord_sf(expand=FALSE) + th +
  labs(title="The covariate z varies WITHIN regions", subtitle="white lines = region borders; peaks sit inside regions")
ggsave("vignettes/t2_raster.png", g1, width=5, height=4.4, dpi=120)

## Fig 2: areal mean choropleth (what naive averaging sees)
g2 <- ggplot(sh) + geom_sf(aes(fill=zbar), color="grey70", linewidth=0.1) +
  scale_fill_viridis_c(name="mean z") + th + labs(title="What areal averaging keeps: the region mean of z")
ggsave("vignettes/t2_zbar.png", g2, width=4.6, height=4, dpi=120)

## Fits: naive GLM vs intensity-scale
naive <- glm(cases ~ zbar + offset(log(pop)), poisson, st_drop_geometry(sh))
fit <- sdalgcp(cases ~ z + offset(log(pop)), data=sh, rasters=r)
sink("vignettes/t2_summary.txt"); print(summary(fit)); sink()

## Fig 3: bias bar chart
bias <- data.frame(method=factor(c("truth","naive areal\naverage","SDALGCP2\n(intensity-scale)"),
                   levels=c("truth","naive areal\naverage","SDALGCP2\n(intensity-scale)")),
                   beta=c(1, coef(naive)["zbar"], fit$beta_opt["z"]))
g3 <- ggplot(bias, aes(method, beta, fill=method)) + geom_col(width=0.6) +
  geom_hline(yintercept=1, linetype="dashed") + geom_text(aes(label=sprintf("%.2f",beta)), vjust=-0.4) +
  scale_fill_manual(values=c("grey60","#B2182B","#2166AC"), guide="none") +
  theme_minimal(base_size=12) + labs(title="Covariate effect: naive vs intensity-scale", x=NULL, y=expression(hat(beta)[z]))
ggsave("vignettes/t2_bias.png", g3, width=5, height=4, dpi=120)

## Fig 4: result maps
pd <- predict(fit, type="discrete")
ggsave("vignettes/t2_arr.png", plot(pd,"adjusted_rr")+th+labs(title="Covariate-adjusted relative risk"), width=4.6,height=4,dpi=120)
cat(sprintf("T2 done: naive=%.2f intensity=%.2f (true 1.0)\n", coef(naive)["zbar"], fit$beta_opt["z"]))
