## Tutorial 5 (spatial confounding / RSR) — a small simulation study + one example.
suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})
make_data <- function() {
  regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(8,8)))
  N <- nrow(regions); cxy <- st_coordinates(st_centroid(regions))
  pts <- sda_points(regions, delta=1.2, method=3)
  S <- as.numeric(t(chol(0.7*precompute_corr(pts, 4)$R[,,1])) %*% rnorm(N))
  regions$x1 <- as.numeric(scale(cxy[,1])); regions$pop <- round(runif(N,500,3000))
  regions$cases <- rpois(N, regions$pop*exp(-6 + 0.6*regions$x1 + S)); regions
}
nrep <- 24; res <- data.frame()
set.seed(100)
for (rr in seq_len(nrep)) {
  d <- make_data()
  std <- sdalgcp(cases~x1+offset(log(pop)), d, control=sdalgcp_control(scale="continuous", n_sim=4000, burnin=1000, thin=5))
  rsr <- sdalgcp(cases~x1+offset(log(pop)), d, control=sdalgcp_control(confounding="restricted", scale="grid", phi=seq(2,7,length.out=5)))
  res <- rbind(res, data.frame(method="standard", est=std$beta_opt["x1"], se=sqrt(std$cov["x1","x1"])),
                    data.frame(method="restricted (RSR)", est=rsr$beta_opt["x1"], se=sqrt(rsr$cov["x1","x1"])))
}
sink("vignettes/t5_study.txt")
agg <- aggregate(cbind(est, se) ~ method, res, function(z) round(c(mean=mean(z), sd=sd(z)),3))
cat("Across", nrep, "simulated data sets (true beta = 0.6):\n")
for (m in unique(res$method)) { sub <- res[res$method==m,]
  cat(sprintf("  %-16s  mean est = %.3f,  SD of est = %.3f,  mean SE = %.3f\n",
      m, mean(sub$est), sd(sub$est), mean(sub$se))) }
sink()
res$method <- factor(res$method, levels=c("standard","restricted (RSR)"))
g <- ggplot(res, aes(method, est, fill=method)) +
  geom_hline(yintercept=0.6, linetype="dashed", color="grey40") +
  geom_boxplot(width=0.5, outlier.size=0.7) +
  scale_fill_manual(values=c("#B2182B","#2166AC"), guide="none") +
  theme_minimal(base_size=12) + labs(title=sprintf("Estimated covariate effect across %d data sets", nrep),
    subtitle="Dashed line = truth (0.6). RSR estimates are tighter and better calibrated.",
    x=NULL, y=expression(hat(beta)[x1]))
ggsave("vignettes/t5_study.png", g, width=5.4, height=4, dpi=120)
cat("T5 done\n")
