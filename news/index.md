# Changelog

## SDALGCP2 0.1.0

First public version.

- **Raster, misaligned, and restricted-regression covariates for the
  spatio-temporal model.**
  [`SDALGCP2_ST()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_ST.md)
  (and `sdalgcp(..., time =)`) now accept `rasters =` (intensity-scale
  continuous covariates), `covariates =` (kriged covariates measured on
  a different support, with the Berkson correction), and
  `confounding = "restricted"` (restricted spatial regression against
  space-time confounding) — the same extensions previously available
  only for spatial fits. The raster/misaligned fits use a Gauss-Newton
  tilting loop around the Kronecker-free space-time likelihood; the
  restricted fit reduces exactly to the spatial restricted fit at
  `T = 1`.

- **Spatio-temporal prediction now returns a long `sf` too.**
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  on an `SDALGCP2_ST` fit returns a one-row-per-region-time `sf` (class
  `"SDALGCP2_ST_pred"`) with the same `relative_risk`/`adjusted_rr`
  columns as the spatial predictor, so both can be mapped or
  [`st_write()`](https://r-spatial.github.io/sf/reference/st_write.html)-en
  the same way; posterior draws are kept as attributes.

- **Bundled datasets.** `sdalgcp_data` — a small simulated `sf` of 64
  regions (`cases`, `x1`, `pop`) used by the help-page examples and the
  intro vignette; and `liver` — a real example, incident primary biliary
  cirrhosis counts by LSOA in North East England (Johnson et al. 2019),
  for realistic case studies.

- **The introductory vignette now runs live on `sdalgcp_data`** (no
  precomputed figures), so it is fully reproducible.

- **Runnable examples on the exported functions**, all using
  `sdalgcp_data`.

- **Prediction output is now an `sf` with clear public-health column
  names.**
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  returns an `sf` (class `"SDALGCP2_pred"`) you can map or
  [`st_write()`](https://r-spatial.github.io/sf/reference/st_write.html)
  directly, with columns `relative_risk`/`relative_risk_se` (the
  relative risk `exp(d'beta + S)`) and `adjusted_rr`/`adjusted_rr_se`
  (the covariate-adjusted relative risk `exp(S)`). The previous
  `RR`/`ARR` names are replaced everywhere
  ([`plot()`](https://rspatial.github.io/terra/reference/plot.html),
  [`exceedance()`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md),
  [`map_exceedance()`](https://olatunjijohnson.github.io/SDALGCP2/reference/map_exceedance.md),
  the spatio-temporal predictor) — `ARR` was dropped because it
  conventionally means *absolute risk reduction* in epidemiology.
  Posterior draws are retained as object attributes so exceedance
  probabilities still work for either quantity.

- **Kronecker-free spatio-temporal model**
  ([`SDALGCP2_ST()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_ST.md)):
  separable space-time SDA-LGCP for counts over the same regions at
  several times. The likelihood never forms the `(N*T)x(N*T)`
  covariance, using `tr(Rs^-1 M Rt^-1 M')` for the quadratic form and
  `N log|Rt| + T log|Rs|` for the log-determinant (both verified vs the
  brute-force Kronecker product). Spatial scale on a grid, temporal
  Matern range estimated continuously (analytic gradient verified vs
  numDeriv). No `geoR`.

- **Covariate-tilted raster model**
  (`SDALGCP2_raster(tilt_spatial = TRUE)`): rebuilds the correlation
  from intensity-tilted weights with a log-normal aggregation
  correction.

- **General Matern smoothness for the direct method**:
  `corr_and_grad_cpp` now provides closed-form phi-derivatives for
  Matern `kappa` in {0.5, 1.5, 2.5}, so continuous-phi
  (`phi_method = "direct"`) estimation works for all three (was
  exponential only). Derivatives verified against finite differences.

- **Nugget / overdispersion term** (`nugget = TRUE`, with
  `phi_method = "direct"`): fits covariance `sigma2 (R(phi) + nu I)` and
  estimates the relative nugget `nu = tau2/sigma2` with a standard
  error, absorbing region-level overdispersion beyond the spatial
  structure. Analytic gradient and Hessian (including the nugget and all
  cross terms) verified against numerical differentiation.

- **Re-anchored (iterated) MCML** (`reanchor =`): re-simulates the
  latent field at the current optimum and refits, keeping the importance
  weights near-uniform. On a 64-region example it lifts the effective
  sample size from ~0% to ~96% and cuts the log-likelihood MC standard
  error ~100x, correcting the variance estimate.

- **Importance-sampling diagnostics**
  ([`mc_diagnostics()`](https://olatunjijohnson.github.io/SDALGCP2/reference/mc_diagnostics.md),
  shown in
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)):
  effective sample size of the importance weights at the optimum and a
  Monte Carlo SE for the maximised log-likelihood; warns on collapse.

- **Spatially continuous (raster) covariates**
  ([`SDALGCP2_raster()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_raster.md)):
  covariates given as rasters enter the LGCP intensity at the
  candidate-point level and are aggregated on the intensity (exp) scale
  via a log-sum-exp offset
  `b_i(beta) = log sum_k w_ik exp(z(x_ik)'beta)` – the correct
  alternative to averaging the predictor over polygons (which is biased
  under the log link). Fit by a Gauss-Newton fixed point reusing
  [`mcml_fit()`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md).
  On a sharp-peak covariate, naive areal averaging is biased +67% while
  this recovers the truth to ~6%.

- **Continuous-phi (“direct”) estimation** (`phi_method = "direct"`):
  optimises the spatial scale `phi` continuously inside the MCML
  objective instead of profiling a grid, using analytic first/second
  derivatives of the aggregated double integral (`corr_and_grad_cpp`).
  Gives a genuine standard error for `phi` from the joint Hessian and
  avoids grid-boundary artefacts; `phi_method = "grid"` stays the
  default. Gradient and Hessian verified against numerical
  differentiation. The full derivation is in
  `math/continuous-phi-derivation.pdf`.

- Post-fit visualisation & diagnostics: relative-risk / uncertainty /
  exceedance maps,
  [`phi_profile()`](https://olatunjijohnson.github.io/SDALGCP2/reference/phi_profile.md),
  [`coef_plot()`](https://olatunjijohnson.github.io/SDALGCP2/reference/coef_plot.md),
  [`model_check()`](https://olatunjijohnson.github.io/SDALGCP2/reference/model_check.md)
  (residual Moran’s I),
  [`report()`](https://olatunjijohnson.github.io/SDALGCP2/reference/report.md).

- C++ (RcppArmadillo + OpenMP) aggregated correlation assembly
  (`precompute_corr`).

- C++ Newton Laplace mode + adaptive MALA sampler (`laplace_sampling`);
  reproducible given the same mode and seed.

- Vectorised, Cholesky-based Monte Carlo likelihood with analytic
  gradient/Hessian (`mcml_fit`).

- One-call spatial fit
  [`SDALGCP2()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2.md)
  (points -\> correlation -\> MCML).

- Candidate-point generation
  [`sda_points()`](https://olatunjijohnson.github.io/SDALGCP2/reference/sda_points.md)
  (SSI / uniform / regular; sf + spatstat

  - terra).

- Prediction
  [`predict.SDALGCP2()`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md)
  (discrete + continuous) with an MCMC or a no-MCMC Laplace fast path;
  [`exceedance()`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md)
  probabilities.

- Vignette and reproducible benchmark scripts using entirely simulated
  data.
