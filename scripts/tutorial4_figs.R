## Tutorial 4 (estimating phi) — grid vs continuous comparison.
suppressMessages({library(SDALGCP2); library(sf); library(ggplot2)})
set.seed(11)
regions <- st_sf(geometry=st_make_grid(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=20,ymax=20))), n=c(9,9)))
N <- nrow(regions); pts <- sda_points(regions, delta=1.1, method=3); cc <- precompute_corr(pts, c(2,3,4,5))
S <- as.numeric(t(chol(0.5*cc$R[,,2])) %*% rnorm(N))            # true phi ~ 3
regions$x1 <- rnorm(N); regions$pop <- round(runif(N,800,4000))
regions$cases <- rpois(N, regions$pop*exp(-6 + 0.5*regions$x1 + S))
phi_grid <- seq(1.5, 6, length.out=12)
tg <- system.time(fg <- sdalgcp(cases~x1+offset(log(pop)), regions,
        control=sdalgcp_control(scale="grid", phi=phi_grid, reanchor=2)))[["elapsed"]]
td <- system.time(fd <- sdalgcp(cases~x1+offset(log(pop)), regions,
        control=sdalgcp_control(scale="continuous", reanchor=2)))[["elapsed"]]
sink("vignettes/t4_compare.txt")
cat(sprintf("GRID:       phi = %.2f            beta_x1 = %.3f   [%.1fs]\n", fg$phi_opt, fg$beta_opt["x1"], tg))
cat(sprintf("CONTINUOUS: phi = %.2f (SE %.2f)  beta_x1 = %.3f   [%.1fs]\n", fd$phi_opt, sqrt(fd$cov["phi","phi"]), fd$beta_opt["x1"], td))
sink()
prof <- data.frame(phi=fg$all_para$phi, dev=-2*(fg$all_para$value-max(fg$all_para$value)))
g <- ggplot(prof, aes(phi, dev)) + geom_line(color="#2166AC", linewidth=1) + geom_point(color="#2166AC") +
  geom_vline(xintercept=fd$phi_opt, color="#B2182B", linewidth=1) +
  annotate("rect", xmin=fd$phi_opt-1.96*sqrt(fd$cov["phi","phi"]), xmax=fd$phi_opt+1.96*sqrt(fd$cov["phi","phi"]),
           ymin=-Inf, ymax=Inf, alpha=0.12, fill="#B2182B") +
  labs(title="Estimating the spatial scale phi", subtitle="Blue: grid profile deviance   Red: continuous estimate (+/-95% CI)",
       x=expression(phi), y="Profile deviance") + theme_minimal(base_size=12)
ggsave("vignettes/t4_compare.png", g, width=6, height=4, dpi=120)
cat("T4 done\n")
