## Tutorial 1 (spatial) — generates figures + captures output. Run from package root.
suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})
th <- theme_minimal(base_size=11) + theme(axis.text=element_blank(), panel.grid=element_blank())
set.seed(2024)
## ---- simulate ----
regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(10,10)))
N <- nrow(regions)
pts <- sda_points(regions, delta=0.9, method=3)
S   <- as.numeric(t(chol(0.5*precompute_corr(pts, 4)$R[,,1])) %*% rnorm(N))
regions$x1    <- rnorm(N)                       # a (non-spatial) covariate
regions$pop   <- round(runif(N, 800, 5000))
regions$cases <- rpois(N, regions$pop*exp(-6 + 0.6*regions$x1 + S))
regions$SIR   <- regions$cases / (regions$pop*exp(-6))   # crude standardised ratio
## ---- data maps ----
gd1 <- ggplot(regions) + geom_sf(aes(fill=SIR), color="grey70", linewidth=0.1) +
  scale_fill_viridis_c(name="crude SIR") + th + labs(title="Observed data: standardised incidence ratio")
gd2 <- ggplot(regions) + geom_sf(aes(fill=x1), color="grey70", linewidth=0.1) +
  scale_fill_distiller(name="x1", palette="RdBu") + th + labs(title="Covariate x1")
ggsave("vignettes/t1_data_sir.png", gd1, width=4.6, height=4, dpi=120)
ggsave("vignettes/t1_data_cov.png", gd2, width=4.6, height=4, dpi=120)
## ---- fit ----
ctrl <- sdalgcp_control(reanchor=3)
fit <- sdalgcp(cases ~ x1 + offset(log(pop)), data=regions, control=ctrl)
sink("vignettes/t1_summary.txt"); print(summary(fit)); sink()
## ---- predictions: discrete RR/ARR/SE/exceedance ----
pd <- predict(fit, type="discrete")
ggsave("vignettes/t1_rr.png",  plot(pd,"relative_risk")  + th + labs(title="Relative risk  exp(eta)"), width=4.6,height=4,dpi=120)
ggsave("vignettes/t1_arr.png", plot(pd,"adjusted_rr") + th + labs(title="Covariate-adjusted RR  exp(S)"), width=4.6,height=4,dpi=120)
ggsave("vignettes/t1_arr_se.png", plot(pd,"adjusted_rr_se") + th + labs(title="SE of adjusted RR"), width=4.6,height=4,dpi=120)
ggsave("vignettes/t1_exc.png", map_exceedance(pd, threshold=1.5, which="adjusted_rr") + th + labs(title="P(adjusted RR > 1.5)"), width=4.6,height=4,dpi=120)
## ---- continuous ----
pc <- predict(fit, type="continuous", sampler="laplace", cellsize=0.6)
ggsave("vignettes/t1_cont.png", plot(pc,"adjusted_rr", bound=regions) + th + labs(title="Continuous adjusted relative risk"), width=4.6,height=4,dpi=120)
## ---- model checking ----
png("vignettes/t1_modelcheck.png", width=520, height=420, res=110)
mc <- model_check(fit, pd); dev.off()
sink("vignettes/t1_check.txt")
cat(sprintf("Residual Moran's I = %.3f (E = %.3f), p = %.3f\n", mc$moran$I, mc$moran$expected, mc$moran$p_value))
d <- mc_diagnostics(fit); cat(sprintf("MC effective sample size: %.0f of %d (%.0f%%)\n", d$ESS, d$B, 100*d$ESS_frac))
sink()
cat("T1 done: x1 =", round(fit$beta_opt["x1"],3), " phi =", round(fit$phi_opt,2), "\n")
