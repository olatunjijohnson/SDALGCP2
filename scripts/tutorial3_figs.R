## Tutorial 3 (spatio-temporal) — figures + output using the new ST plot method.
suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})
th <- theme_minimal(base_size=11) + theme(axis.text=element_blank(), panel.grid=element_blank())
set.seed(7)
shp <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=18,ymax=18))), n=c(6,6)))
N <- nrow(shp); times <- 2019:2022; T <- length(times)
pts <- sda_points(shp, delta=1.4, method=3); cc <- precompute_corr(pts, c(2,3,4))
Rs <- cc$R[,,2]; Rt <- SDALGCP2:::.temporal_corr(seq_len(T), 1.5, 0.5)
L <- t(chol(0.4*kronecker(Rt,Rs)))
x1 <- rnorm(N*T); pop <- round(runif(N*T,1000,5000))
y <- rpois(N*T, pop*exp(-6 + 0.5*x1 + as.numeric(L%*%rnorm(N*T))))
dat <- st_sf(data.frame(cases=y, x1=x1, pop=pop, year=rep(times,each=N)),
             geometry=st_geometry(shp)[rep(seq_len(N),T)])
## data map: SIR by year
dat$SIR <- dat$cases/(dat$pop*exp(-6))
gd <- ggplot(dat) + geom_sf(aes(fill=SIR), color="grey70", linewidth=0.1) + facet_wrap(~year, nrow=1) +
  scale_fill_viridis_c(name="crude\nSIR") + th + labs(title="Observed data by year (crude SIR)")
ggsave("vignettes/t3_data.png", gd, width=8, height=2.6, dpi=120)
## fit
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data=dat, time="year",
               control=sdalgcp_control(reanchor=2, n_sim=6000, burnin=1500, thin=5))
sink("vignettes/t3_summary.txt"); print(summary(fit)); sink()
## predict + maps with the new API
pr <- predict(fit)
ggsave("vignettes/t3_arr_facet.png", plot(pr, time=NULL, what="adjusted_rr")+th+labs(title="Covariate-adjusted relative risk by year"),
       width=8, height=2.6, dpi=120)
ggsave("vignettes/t3_rr_2021.png", plot(pr, time=2021, what="relative_risk")+th, width=4.4, height=4, dpi=120)
ggsave("vignettes/t3_exc_2021.png", plot(pr, time=2021, what="exceedance", threshold=1.3, which="adjusted_rr")+th, width=4.4, height=4, dpi=120)
cat(sprintf("T3 done: x1=%.3f phi=%.2f nu=%.2f\n", fit$beta_opt["x1"], fit$phi_opt, fit$nu_opt))
