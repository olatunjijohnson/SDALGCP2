suppressMessages({library(SDALGCP2); library(sf); library(terra); library(ggplot2)})
th <- theme_minimal(base_size=11) + theme(axis.text=element_blank(), panel.grid=element_blank())

## ---------- Spatio-temporal facet map ----------
set.seed(7)
shp <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=18,ymax=18))), n=c(6,6)))
N <- nrow(shp); T <- 4; times <- 2019:2022
pts <- sda_points(shp, delta=1.4, method=3); cc <- precompute_corr(pts, c(2,3,4))
Rs <- cc$R[,,2]; Rt <- SDALGCP2:::.temporal_corr(seq_len(T), 1.5, 0.5)
L <- t(chol(0.4*kronecker(Rt,Rs)))
x1 <- rnorm(N*T); pop <- round(runif(N*T,1000,5000))
y <- rpois(N*T, pop*exp(-6 + 0.5*x1 + as.numeric(L%*%rnorm(N*T))))
dat <- st_sf(data.frame(cases=y, x1=x1, pop=pop, year=rep(times,each=N)),
             geometry=st_geometry(shp)[rep(seq_len(N),T)])
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data=dat, time="year",
               control=sdalgcp_control(reanchor=2, n_sim=5000, burnin=1500, thin=5))
pred <- predict(fit)
map_sf <- do.call(rbind, lapply(seq_len(T), function(t)
  { g <- shp; g$ARR <- pred$ARR_mean[,t]; g$year <- times[t]; g }))
gst <- ggplot(map_sf) + geom_sf(aes(fill=ARR), color="grey70", linewidth=0.1) +
  facet_wrap(~year, nrow=1) +
  scale_fill_gradient2(name="Relative\nrisk", midpoint=1, low="#2166AC", mid="grey95", high="#B2182B") +
  th + labs(title="Spatio-temporal relative risk by year")
ggsave("vignettes/figures/st_facet.png", gst, width=8, height=2.6, dpi=120)
cat("ST done: beta x1 =", round(fit$beta_opt["x1"],3), "\n")

## ---------- Raster: naive vs intensity-scale bias ----------
set.seed(3)
r <- rast(xmin=0,xmax=20,ymin=0,ymax=20,resolution=0.1); xy <- xyFromCell(r,1:ncell(r))
src <- matrix(runif(2*14,1,19),ncol=2)
values(r) <- rowSums(sapply(1:nrow(src),function(s)3.2*exp(-((xy[,1]-src[s,1])^2+(xy[,2]-src[s,2])^2)/0.5)))
names(r)<-"z"
sh <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))),n=c(8,8)))
Nr<-nrow(sh); pp<-sda_points(sh,delta=0.5,method=3)
Z<-lapply(pp,function(p)cbind(1,terra::extract(r,as.matrix(p$xy))[,"z"])); w<-lapply(pp,function(p)p$weight)
bt<-c(-6,1); bi<-sapply(1:Nr,function(i)log(sum(w[[i]]*exp(as.numeric(Z[[i]]%*%bt)))))
popr<-round(runif(Nr,2000,6000)); yr<-rpois(Nr,popr*exp(bi))
zbar<-sapply(1:Nr,function(i)sum(w[[i]]*Z[[i]][,2]))
naive<-glm(yr~zbar+offset(log(popr)),poisson)
fr <- sdalgcp(cases ~ z + offset(log(pop)),
              data=st_sf(data.frame(cases=yr,pop=popr),geometry=st_geometry(sh)), rasters=r)
bias <- data.frame(method=c("truth","naive areal\naverage","SDALGCP2\n(intensity-scale)"),
                   beta=c(1, coef(naive)["zbar"], fr$beta_opt["z"]))
bias$method <- factor(bias$method, levels=bias$method)
gb <- ggplot(bias, aes(method, beta, fill=method)) + geom_col(width=0.6) +
  geom_hline(yintercept=1, linetype="dashed") +
  geom_text(aes(label=sprintf("%.2f",beta)), vjust=-0.4) +
  scale_fill_manual(values=c("grey60","#B2182B","#2166AC"), guide="none") +
  theme_minimal(base_size=12) + labs(title="Covariate effect: naive vs intensity-scale", x=NULL, y=expression(hat(beta)[z]))
ggsave("vignettes/figures/raster_bias.png", gb, width=5, height=4, dpi=120)
cat("raster done: naive=", round(coef(naive)["zbar"],2), " intensity=", round(fr$beta_opt["z"],2), "\n")
