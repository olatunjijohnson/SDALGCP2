suppressMessages({library(SDALGCP2); library(SDALGCP)})
set.seed(2024)
N <- 60                      # regions
centers <- matrix(runif(N*2,0,30), ncol=2)
points <- lapply(seq_len(N), function(i){
  ni <- sample(6:12,1)
  xy <- sweep(matrix(rnorm(ni*2,0,1.2),ncol=2), 2, centers[i,], "+")
  w  <- rep(1/ni, ni)
  list(xy=xy, weight=w)
})
attr(points,"weighted") <- FALSE
attr(points,"my_shp")  <- NULL
phi_grid <- seq(1.5, 6, length.out = 10)
corr <- precompute_corr(points, phi_grid)              # shared correlation array

## simulate areal counts from the SDA model at known params
phi_t <- 3.0; sig2_t <- 0.6; beta_t <- c(-0.2, 0.8)
k <- which.min(abs(phi_grid - phi_t))
Sig <- sig2_t * corr$R[,,k]
x1 <- rnorm(N)
D  <- cbind(1, x1)
Strue <- as.numeric(t(chol(Sig)) %*% rnorm(N))
m  <- round(runif(N, 50, 200))
y  <- rpois(N, m*exp(D %*% beta_t + Strue))
dat <- data.frame(y=y, x1=x1, pop=m)
FORM <- y ~ x1 + offset(log(pop))
ctrl <- list(n.sim=8000, burnin=2000, thin=6, h=1.65/N^(1/6), c1.h=0.01, c2.h=1e-4)

set.seed(99)
t2 <- system.time(fit2 <- mcml_fit(FORM, dat, corr, control.mcmc=ctrl))[["elapsed"]]
set.seed(99)
t1 <- system.time(fit1 <- SDALGCP:::SDALGCPParaEst(FORM, dat, corr, control.mcmc=ctrl))[["elapsed"]]

cat("True:    beta=",paste(round(beta_t,3),collapse=", "),
    " sigma2=",sig2_t," phi=",phi_t,"\n",sep="")
cat("SDALGCP2 beta=",paste(round(fit2$beta_opt,3),collapse=", "),
    " sigma2=",round(fit2$sigma2_opt,3)," phi=",round(fit2$phi_opt,3),
    sprintf("  (%.2fs)\n",t2),sep="")
cat("SDALGCP  beta=",paste(round(as.numeric(fit1$beta_opt),3),collapse=", "),
    " sigma2=",round(fit1$sigma2_opt,3)," phi=",round(fit1$phi_opt,3),
    sprintf("  (%.2fs)\n",t1),sep="")
cat(sprintf("Speedup (estimation): %.1fx\n", t1/t2))
cat("\n--- S3 methods ---\n")
print(summary(fit2))
cat("\nconfint:\n"); print(confint(fit2))
