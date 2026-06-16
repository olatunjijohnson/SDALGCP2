# SDALGCP2 0.0.0.9000 (development)

First development version — a faster, modernised successor to SDALGCP.

* **Nugget / overdispersion term** (`nugget = TRUE`, with `phi_method = "direct"`):
  fits covariance `sigma2 (R(phi) + nu I)` and estimates the relative nugget
  `nu = tau2/sigma2` with a standard error, absorbing region-level overdispersion
  beyond the spatial structure. Analytic gradient and Hessian (including the nugget
  and all cross terms) verified against numerical differentiation.
* **Re-anchored (iterated) MCML** (`reanchor =`): re-simulates the latent field at
  the current optimum and refits, keeping the importance weights near-uniform. On a
  64-region example it lifts the effective sample size from ~0% to ~96% and cuts the
  log-likelihood MC standard error ~100x, correcting the variance estimate.
* **Importance-sampling diagnostics** (`mc_diagnostics()`, shown in `summary()`):
  effective sample size of the importance weights at the optimum and a Monte Carlo
  SE for the maximised log-likelihood; warns on collapse.
* **Spatially continuous (raster) covariates** (`SDALGCP2_raster()`): covariates
  given as rasters enter the LGCP intensity at the candidate-point level and are
  aggregated on the intensity (exp) scale via a log-sum-exp offset
  `b_i(beta) = log sum_k w_ik exp(z(x_ik)'beta)` -- the correct alternative to
  averaging the predictor over polygons (which is biased under the log link). Fit
  by a Gauss-Newton fixed point reusing `mcml_fit()`. On a sharp-peak covariate,
  naive areal averaging is biased +67% while this recovers the truth to ~6%.
* **Continuous-phi ("direct") estimation** (`phi_method = "direct"`): optimises the
  spatial scale `phi` continuously inside the MCML objective instead of profiling a
  grid, using analytic first/second derivatives of the aggregated double integral
  (`corr_and_grad_cpp`). Gives a genuine standard error for `phi` from the joint
  Hessian and avoids grid-boundary artefacts; `phi_method = "grid"` stays the
  default. Gradient and Hessian verified against numerical differentiation. The
  full derivation is in `math/continuous-phi-derivation.pdf`.
* Post-fit visualisation & diagnostics: relative-risk / uncertainty / exceedance
  maps, `phi_profile()`, `coef_plot()`, `model_check()` (residual Moran's I),
  `report()`.
* C++ (RcppArmadillo + OpenMP) aggregated correlation assembly (`precompute_corr`).
* C++ Newton Laplace mode + adaptive MALA sampler (`laplace_sampling`); bit-identical
  to `SDALGCP::Laplace.sampling` given the same mode and seed.
* Vectorised, Cholesky-based Monte Carlo likelihood with analytic gradient/Hessian
  (`mcml_fit`); matches the original to ~1e-14 and is ~10x faster on estimation.
* One-call spatial fit `SDALGCP2()` (points -> correlation -> MCML); reproduces
  `SDALGCP::SDALGCPMCML` estimates and is ~8x faster end-to-end on a 64-region example.
* Candidate-point generation `sda_points()` (SSI / uniform / regular; sf + spatstat
  + terra; drops splancs/sp).
* Prediction `predict.SDALGCP2()` (discrete + continuous) with an MCMC or a no-MCMC
  Laplace fast path; `exceedance()` probabilities.
* Vignette and reproducible benchmark scripts using entirely simulated data.
