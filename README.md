# SDALGCP2

A faster, modernised successor to
[**SDALGCP**](https://github.com/olatunjijohnson/SDALGCP) for fitting a
**S**patially **D**iscrete **A**pproximation to a **L**og-**G**aussian **C**ox
**P**rocess (SDA-LGCP) to spatially aggregated disease counts.

SDALGCP2 keeps the statistical method of Johnson, Diggle & Giorgi (2019) but

- moves the performance-critical kernels (aggregated correlation assembly, MALA
  sampling, the Monte Carlo likelihood, and the spatio-temporal Kronecker
  likelihood) into **C++ via RcppArmadillo + OpenMP**;
- removes orphaned/legacy dependencies (`geoR`, `sp`, `raster`, `spacetime`,
  `mapview`, `splancs`, `pdist`, `maxLik`) in favour of `sf`/`terra`/`stars`;
- adds statistical options: a **nugget** term, general **Matern** smoothness,
  importance-sampling **diagnostics** (effective sample size, MC error), and
  **re-anchored** MCML for robustness.

See [`DESIGN.md`](DESIGN.md) for the full analysis of the original package, the
bottlenecks, and the build plan.

> **Status:** scaffolding + flagship correlation kernel implemented. This is a
> work in progress; see the build order in `DESIGN.md`.

## Install (development)

```r
# from the parent directory containing SDALGCP2/
devtools::install("SDALGCP2")
```

## What works today

```r
library(SDALGCP2)
# `points` is the candidate-point list from the point-generation step
corr <- precompute_corr(points, phi = seq(500, 1700, length.out = 20))
```
