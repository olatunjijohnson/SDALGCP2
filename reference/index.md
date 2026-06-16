# Package index

## Model fitting

User entry points for fitting SDA-LGCP models.

- [`SDALGCP2()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2.md)
  : Fit a spatial SDA-LGCP model
- [`SDALGCP2_raster()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_raster.md)
  : Fit an SDA-LGCP with spatially continuous (raster) covariates
- [`SDALGCP2_ST()`](https://olatunjijohnson.github.io/SDALGCP2/reference/SDALGCP2_ST.md)
  : Fit a spatio-temporal SDA-LGCP model (Kronecker-free)
- [`mcml_fit()`](https://olatunjijohnson.github.io/SDALGCP2/reference/mcml_fit.md)
  : Monte Carlo maximum likelihood estimation for the spatial SDA-LGCP

## Candidate points and correlation

Discretisation and the aggregated correlation assembly.

- [`sda_points()`](https://olatunjijohnson.github.io/SDALGCP2/reference/sda_points.md)
  : Generate candidate sampling points inside each region
- [`precompute_corr()`](https://olatunjijohnson.github.io/SDALGCP2/reference/precompute_corr.md)
  : Precompute aggregated region-level correlation matrices

## Sampling and control

- [`control_mcmc()`](https://olatunjijohnson.github.io/SDALGCP2/reference/control_mcmc.md)
  : MCMC control settings for the MALA sampler
- [`laplace_sampling()`](https://olatunjijohnson.github.io/SDALGCP2/reference/laplace_sampling.md)
  : Sample the latent field \[S \| Y\] (Poisson, non-nested) via C++
  MALA

## Prediction

- [`predict(`*`<SDALGCP2>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2.md)
  : Predict relative risk from a fitted SDALGCP2 model
- [`predict(`*`<SDALGCP2_ST>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/predict.SDALGCP2_ST.md)
  : Discrete (region x time) prediction for a spatio-temporal fit
- [`exceedance()`](https://olatunjijohnson.github.io/SDALGCP2/reference/exceedance.md)
  : Exceedance probabilities P(relative risk \> threshold)

## Visualisation

Maps and post-fit graphics.

- [`plot(`*`<SDALGCP2_pred>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/plot.SDALGCP2_pred.md)
  : Map a fitted SDALGCP2 prediction
- [`plot(`*`<SDALGCP2>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/plot.SDALGCP2.md)
  : Plot an SDALGCP2 fit (the phi profile deviance)
- [`map_exceedance()`](https://olatunjijohnson.github.io/SDALGCP2/reference/map_exceedance.md)
  : Map exceedance probabilities P(relative risk \> threshold)
- [`phi_profile()`](https://olatunjijohnson.github.io/SDALGCP2/reference/phi_profile.md)
  : Profile likelihood and confidence interval for the spatial scale phi
- [`coef_plot()`](https://olatunjijohnson.github.io/SDALGCP2/reference/coef_plot.md)
  : Coefficient plot of fixed effects (and sigma^2) with confidence
  intervals
- [`report()`](https://olatunjijohnson.github.io/SDALGCP2/reference/report.md)
  : One-call panel of post-fit graphics

## Diagnostics

- [`mc_diagnostics()`](https://olatunjijohnson.github.io/SDALGCP2/reference/mc_diagnostics.md)
  : Importance-sampling diagnostics for an MCML fit
- [`model_check()`](https://olatunjijohnson.github.io/SDALGCP2/reference/model_check.md)
  : Posterior-predictive model checking for an SDALGCP2 fit

## Methods

- [`summary(`*`<SDALGCP2>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/summary.SDALGCP2.md)
  : Summary of an SDALGCP2 fit
- [`confint(`*`<SDALGCP2>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/confint.SDALGCP2.md)
  : Wald confidence intervals for an SDALGCP2 fit
- [`print(`*`<SDALGCP2>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/print.SDALGCP2.md)
  : Print an SDALGCP2 fit
- [`print(`*`<summary.SDALGCP2>`*`)`](https://olatunjijohnson.github.io/SDALGCP2/reference/print.summary.SDALGCP2.md)
  : Print a summary of an SDALGCP2 fit
