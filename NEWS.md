# SDALGCP2 0.0.0.9000 (development)

First development version — a faster, modernised successor to SDALGCP.

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
