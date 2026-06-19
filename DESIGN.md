# SDALGCP2 — Design Notes

A fast, statistically rich implementation of the spatially discrete approximation
to a log-Gaussian Cox process.

> SDALGCP2 fits a *spatially discrete approximation* to a log-Gaussian Cox process
> (SDA-LGCP) for spatially aggregated disease counts, estimating parameters by
> Monte Carlo Maximum Likelihood (MCML) and predicting discrete and continuous
> relative risk (method: Johnson, Diggle & Giorgi, 2019, *Stat. Med.* 38:4871–4887).

This document records (1) what the method does, (2) where the current package is
slow or statistically weak, and (3) the concrete plan for SDALGCP2.

---

## 1. The model and the approximation (so the speed work stays faithful)

Disease counts `Y_i` in regions `A_i`, `i = 1..N`, are modelled as a
log-Gaussian Cox process aggregated to area level:

```
Y_i | S  ~  Poisson( m_i * exp( d_i' beta + T_i(S) ) )
S(.)      ~  zero-mean Gaussian process,  Cov{S(x),S(x')} = sigma2 * exp(-||x-x'||/phi)
```

The continuous areal integral `∫_{A_i} exp(S(x)) dx` is replaced by a **discrete
average over candidate points** placed inside each region:

```
T_i(S) = sum_k w_ik S(x_ik)         (population-weighted, sum_k w_ik = 1)
       = mean_k S(x_ik)             (unweighted)
```

Hence the *region-level* random effect vector `S* = (T_1,...,T_N)` is Gaussian
with an `N x N` covariance

```
Sigma(phi, sigma2)_{ij} = sigma2 * sum_{k,l} w_ik w_jl exp(-||x_ik - x_jl|| / phi)
```

Everything downstream (MCML, prediction) operates on this `N x N` object. The
candidate-point machinery only exists to *build* `Sigma`. **This split — an
expensive one-off "build the aggregated correlation" step, then repeated cheap
`N x N` linear algebra — is the key to where time is spent.**

Estimation (Christensen 2004 MCML):
1. Pick an anchor `(phi0, sigma2_0)`, draw `S*` samples from `[S* | Y]` by a
   Langevin–Hastings (MALA) sampler.
2. For each `phi` on a grid, maximise an importance-sampling-corrected Monte
   Carlo log-likelihood over `(beta, sigma2)` with `phi` fixed.
3. Profile over `phi`; report the maximiser and a deviance-based CI.

Spatio-temporal version: separable `Sigma_space ⊗ Sigma_time`, time correlation
Matérn(`kappa`) with range `nu`.

---

## 2. Where a naive implementation is slow (measured against the algorithm above)

Ranked by impact.

### B1. `precomputeCorrMatrix` — dominant cost
Pure-R double loop over all `N(N+1)/2` region pairs; for each pair it forms an
`n_i x n_j` distance matrix and a **3-D array** `exp(-outer(D, 1/phi))`, then
reduces it for every `phi`. Cost ≈ `O(N^2 * nbar^2 * n_phi)` with large constant
factors from R-level `outer`, array allocation, and `pdist`.
- For the PBC example (N = 545, ~tens of points/region, 20 `phi`) this is the
  wall-clock bottleneck and allocates gigabytes transiently.
- **Fix:** single RcppArmadillo kernel, OpenMP-parallel over region pairs.
  Compute each pairwise distance **once**, accumulate the weighted/mean
  exponential sum for all `phi` in a tight inner loop — no 3-D array, no per-phi
  R allocation. Expected 1–2 orders of magnitude faster.

### B2. `Laplace.sampling` — the MALA loop is interpreted R
`n.sim` (default 10 000) iterations, each doing `Sigma.sroot %*% W` and gradient
products in R. Runs once for estimation and again for *every* prediction call.
- **Fix:** port the inner loop to C++ (the proposal density, gradient, accept
  step). Keep adaptation identical so results match. ~10–30x on the sampler.

### B3. `Aggregated_poisson_log_MCML` — un-vectorised MC likelihood
`apply(S.sim, 1, ...)` plus `lapply/Reduce` for gradient and Hessian rebuild the
same quadratic forms sample-by-sample on every optimiser step.
- **Fix:** vectorise. With `S.sim` (`n.sim x N`) and a Cholesky factor of `R`,
  all quadratic forms `(S-mu)' R^{-1} (S-mu)` are one `crossprod`/triangular
  solve; gradients/Hessians become matrix ops. Factor `R` once (Cholesky →
  log-det + solves), never call `solve()` in the inner loop.

### B4. Spatio-temporal Kronecker blow-up
`Aggregated_poisson_log_MCML_ST` materialises `kronecker(inv.t.corr, inv.s.corr)`
— an `(N*T) x (N*T)` dense matrix — and its log-det/products **on every
likelihood, gradient and Hessian call**. With Kronecker structure this is
unnecessary:
```
(A ⊗ B)^{-1} = A^{-1} ⊗ B^{-1}
det(A ⊗ B)   = det(A)^{n_B} det(B)^{n_A}
x'(A⊗B)^{-1}x = tr( B^{-1} X' A^{-1} X ),   X = mat(x)  (T x N or N x T)
```
- **Fix:** never form the Kronecker product; reshape each sample to a `T x N`
  matrix and use the identities. Turns `(N*T)^3` work into `N^3 + T^3` plus
  cheap matrix multiplies. This is the single biggest ST win.
- Also drop `geoR::varcov.spatial` (orphaned dependency) for the temporal Matérn
  — compute Matérn directly in a few lines / C++.

### B5. `compute_cross_cov` (continuous prediction) — R double loop
`n.pred * N` iterations each calling `pdist`. **Fix:** C++ kernel mirroring B1.

### B6. Pervasive `solve()` / `determinant()` instead of Cholesky
Many sites invert full matrices and separately compute determinants. Replace
with one Cholesky per matrix → triangular solves + `2*sum(log(diag(L)))`.

### B7. `maxim.integrand` uses generic `maxBFGS`
Analytic gradient **and** Hessian are available and the Hessian is
`-Sigma^{-1} - diag(h)`. A handful of damped Newton steps converge faster and
more reliably than BFGS. **Fix:** custom Newton with the supplied Hessian.

---

## 3. Statistical improvements (thinking as a statistician)

### S1. General Matérn spatial covariance (not just exponential)
The point-level kernel is hard-wired to `exp(-d/phi)` (Matérn `kappa=0.5`).
Expose `kappa` (½, 1½, 2½ closed forms; general via Bessel) so smoothness is a
modelling choice, consistent with the temporal side which already uses Matérn.

### S2. A nugget / unstructured area effect (overdispersion)
Add an optional region-level i.i.d. term `U_i ~ N(0, tau2)` so
`Sigma = sigma2 * R(phi) + tau2 * I`. Aggregated disease counts routinely show
overdispersion beyond the spatial signal; without a nugget `sigma2`/`phi` absorb
it and bias the range. This is a small change to the covariance assembly and a
genuinely useful modelling option.

### S3. Honest Monte Carlo error & importance-sampling diagnostics
The MCML estimate is conditional on draws taken at the anchor `phi0`. Currently
there is **no** diagnostic for how trustworthy the importance reweighting is when
`phi` is far from `phi0`.
- Report the **effective sample size** of the importance weights at the chosen
  `phi`, `ESS = (sum w)^2 / sum w^2`, and warn when it collapses.
- Report a Monte Carlo standard error for the maximised log-likelihood.

### S4. Iterated / re-anchored MCML
After the first profile, **re-anchor** the sampler at the current `(phi, sigma2)`
optimum and refit once or twice (stop when estimates stabilise). This fixes the
main failure mode of single-pass MCML — a wide `phi` grid far from the truth —
at modest extra cost, and makes results far less sensitive to `par0`.

### S5. Proper profile CI for `phi`
Replacing the MC log-likelihood profile with a **loess** smooth before taking the
deviance CI (current `phiCI`) smooths Monte Carlo noise *and* genuine curvature
together. Prefer a monotone spline through the (MC-error-aware) profile points,
or a finer adaptive `phi` grid near the maximum, and propagate MC error into the
reported interval.

### S6. Optional Laplace (no-MCMC) fast path for prediction
For Poisson-aggregated data the Gaussian/Laplace approximation to `[S*|Y]` is
typically excellent. Offer `method = "laplace"` that skips MALA entirely for
prediction (mode + Gaussian covariance already computed in `maxim.integrand`).
Near-instant maps; keep MCMC as the accurate default.

### S7. Parameter-uncertainty propagation in prediction
Continuous/discrete prediction plugs in MCML point estimates. Offer an option to
propagate `(beta, sigma2, phi)` uncertainty (draw from the asymptotic / profile
distribution) so predictive SEs are not over-confident.

### S8. Data-driven default `phi` grid
Current default `[sqrt(min area), extent/10]` is geometric and arbitrary. Seed
the grid from an empirical variogram of Poisson-GLM Pearson residuals so the grid
brackets a sensible range; still allow manual override.

---

## 3b. Methodology: spatially continuous (raster) predictors

**Goal (user request):** let covariates be supplied as rasters `z(x)` and enter the
model *correctly*, **not** by first averaging the predictor over each polygon and
plugging the average into a region-level GLM. Naive areal aggregation
`(\bar z_i)' \beta` is biased whenever the exp link is nonlinear and `z` varies
within a region — it ignores Jensen's inequality.

**Principle:** aggregate on the *intensity* scale, not the linear-predictor scale.
The aggregated-LGCP mean for region `i` is
```
E[Y_i | S] = m_i ∫_{A_i} exp(z(x)'beta + S(x)) rho(x) dx
           ≈ m_i * sum_k w_ik exp(z(x_ik)'beta + S(x_ik))          (discretised)
```
where `x_ik` are the candidate points already used to build `Sigma` and `w_ik` are
the (population) weights. Factor the covariate part out:
```
= m_i * exp(b_i(beta)) * sum_k c_ik(beta) exp(S(x_ik)),
  b_i(beta)  = log sum_k w_ik exp(z(x_ik)'beta)        # log-sum-exp covariate offset
  c_ik(beta) = w_ik exp(z(x_ik)'beta) / exp(b_i(beta)) # covariate-tilted weights, sum_k = 1
```
Then approximate the inner sum by its log-normal aggregate (linear functional of `S`
plus a variance correction):
```
sum_k c_ik exp(S(x_ik)) ≈ exp( T_i(S;beta) + 0.5 v_i(beta) ),
  T_i(S;beta) = sum_k c_ik(beta) S(x_ik)               # covariate-tilted spatial avg
  v_i(beta)   = Var of that aggregation (from the point-level covariance)
```
giving a model that **keeps the fast N-dimensional structure**:
```
E[Y_i | S*] = m_i * exp( b_i(beta) + T_i(S*;beta) + 0.5 v_i(beta) ),
S* ~ N(0, Sigma(phi, sigma2; beta)),  Sigma built with the tilted weights c_ik.
```

**Why this is the right contribution**
- `b_i(beta)` is the **log-sum-exp** of the within-region covariate surface — it is
  *not* `(\bar z_i)'beta`, and reduces to it only when `z` is constant in `A_i`.
- The spatial aggregation weights `c_ik(beta)` are themselves tilted by the
  covariate intensity, so high-covariate sub-areas contribute more to the region
  effect — the statistically coherent behaviour.
- It collapses to the basic region-constant model exactly when covariates are
  region-constant (`c_ik = w_ik`, `b_i = z_i'beta`, `v_i` absorbed into the intercept).

**Computational consequence & plan.** `Sigma` and the offset now depend on `beta`,
so the clean "precompute correlation once, profile cheaply" loop becomes an
*iterated* scheme: (i) hold `beta` fixed → rebuild tilted weights → assemble
`Sigma(beta)` and `b(beta)` via the existing C++ kernel (which already takes
arbitrary weights), (ii) run the MCML step, (iii) update `beta`, iterate to
convergence. A cheaper **first-order** option keeps `c_ik = w_ik` (covariate enters
only through the offset `b_i(beta)` and `v_i`), leaving `Sigma` `beta`-independent
and fully precomputable — still strictly better than areal averaging. Both will be
offered; default to the first-order offset model, with the fully-tilted version as
an option. Implementation tracked as task #5; design fixed here.

## 3c. Methodology: continuous-phi (grid-free) MCML — "solve it directly"

**Current behaviour.** The spatial scale `phi` is *discretised*: we evaluate the
aggregated correlation `R(phi)` on a grid, fit `(beta, sigma2)` at each, and take
the profile maximum. The grid is robust (it shows the whole profile) but wasteful
— most grid points are far from the optimum — and the resolution caps the
precision of `phi-hat`.

**Idea (user).** Treat `phi` as a continuous parameter and optimise it *directly*
inside the MCML objective, alongside `beta` and `sigma2`, instead of scanning a
grid. The MCML log-likelihood is smooth in `phi`, so a gradient-based optimiser
needs only `R(phi)` and its derivative `dR/dphi` at arbitrary `phi`.

**Derivatives are available in closed form.** For the aggregated exponential
kernel,
```
R_ij(phi)      = sum_{k,l} w_ik w_jl exp(-d_kl/phi)
dR_ij/dphi     = sum_{k,l} w_ik w_jl exp(-d_kl/phi) * (d_kl / phi^2)
```
(general Matern: differentiate the closed-form 3/2, 5/2 kernels similarly). The
MCML gradient w.r.t. `log phi` then follows exactly as in the existing
spatio-temporal code, which already does this for the temporal range `nu`:
```
d/dphi of the per-sample joint log-density
  = -0.5 [ tr(R^{-1} R')  -  (S-mu)' R^{-1} R' R^{-1} (S-mu) / sigma2 ],   R' = dR/dphi
```
weighted by the importance weights and reduced over samples — the same vectorised
machinery as `mcml_fit`, plus one extra parameter.

**Plan.**
- New C++ kernel `corr_and_grad_cpp(coords, weights, phi, kappa)` returning `R`
  and `dR/dphi` (reuses the `corr_aggregate` inner loop).
- `mcml_fit(..., method = c("grid","direct"))`: `"grid"` stays the default and the
  reference; `"direct"` optimises `(beta, log sigma2, log phi)` jointly with
  `nlminb`/`optim`, recomputing `R(phi)`/`R'(phi)` per step (cheap now).
- Standard errors for `phi` come directly from the joint Hessian (no loess/spline
  needed), and the profile plot can still be produced on demand.
- **Compare** the two: agreement of `phi-hat` and speed (direct should need far
  fewer correlation builds than a 20-point grid). Keep both options permanently —
  grid for diagnostics/multimodal profiles, direct for speed/precision.

Tracked as task #6; design fixed here.

## 4. Engineering / dependency modernisation

- **Drop legacy/orphaned deps:** `geoR` (orphaned), `sp`, `raster`, `spacetime`,
  `mapview`, `splancs`, `pdist`, `maxLik`. Keep a lean stack: `sf`, `terra`,
  `stars`, `spatstat.geom`/`spatstat.random` (sampling), `Matrix`, `ggplot2`,
  plus `Rcpp`/`RcppArmadillo` for the kernels.
- **Pitfalls to avoid:**
  - `class(para_est) != "..."` — exact class comparison breaks with multi-class
    objects; use `inherits()`.
  - `scale_fill_viridis_c(name = cat(...))` — `cat()` returns `NULL`; the legend
    title is silently dropped.
  - `aes_string()` is deprecated — use tidy-eval `.data[[var]]`.
  - `SDADiscretePred` recomputes `apply(.,2,sd)` repeatedly on `exp(S.sim)`; cache.
- **Reproducibility:** thread an explicit RNG seed through samplers; document that
  OpenMP reductions are order-stable for the correlation build.
- **Testing:** `testthat` suite asserting SDALGCP2 reproduces SDALGCP estimates on
  the PBC data to Monte-Carlo tolerance (correctness gate before/while optimising).
- **Benchmarks:** `bench`-based scripts comparing each kernel old vs new.

---

## 5. Proposed package architecture

```
SDALGCP2/
  DESCRIPTION            # lean deps + LinkingTo Rcpp, RcppArmadillo
  NAMESPACE              # roxygen-generated
  R/
    covariance.R         # aggregated correlation assembly (R wrapper over C++), Matern, nugget
    sampling.R           # Laplace mode (Newton) + MALA (C++-backed)
    mcml.R               # vectorised MC likelihood, profile over phi, re-anchoring
    fit.R                # SDALGCP2() user entry point (spatial)
    fit_st.R             # spatio-temporal (Kronecker-free)
    predict.R            # discrete + continuous prediction, Laplace fast path
    points.R             # candidate-point generators (SSI / uniform / regular)
    methods.R            # print/summary/confint/plot, exceedance
    diagnostics.R        # ESS, MC error, profile CI for phi
  src/
    corr_aggregate.cpp   # B1: OpenMP aggregated correlation (FLAGSHIP, included)
    mala_sampler.cpp     # B2
    mc_loglik.cpp        # B3
    kron_st.cpp          # B4 helpers
    cross_cov.cpp        # B5
  tests/testthat/
  vignettes/
  DESIGN.md              # this file
```

### Naming
Keep user-facing verbs recognisable and consistent:
`SDALGCP2()` (fit), `predict()` S3 method, `control_mcmc()`, `summary()`,
`confint()`, `phi_profile()`, `plot()`.

---

## 6. Expected impact (rough, to be confirmed by benchmarks)

| Stage                         | Current        | SDALGCP2 target            |
|-------------------------------|----------------|----------------------------|
| Correlation precompute (B1)   | minutes        | seconds (10–50x, parallel) |
| MALA sampler (B2)             | tens of sec    | ~10–30x                    |
| MC likelihood profile (B3)    | seconds–min    | ~5–20x (vectorised)        |
| ST likelihood (B4)            | infeasible-ish | `N^3+T^3` not `(NT)^3`     |
| Continuous prediction (B5)    | minutes        | seconds                    |

Numbers are targets; the test suite gates correctness first, then `bench`
confirms the speedups.

---

## 7. Build order

1. **Scaffold + flagship kernel** (this commit): package skeleton, DESIGN.md,
   `corr_aggregate.cpp` + R wrapper + a correctness/benchmark check vs the
   reference R implementation.
2. Vectorised MC likelihood (B3) + Newton Laplace mode (B7).
3. C++ MALA (B2), wire spatial fit end-to-end, `testthat` correctness checks.
4. Prediction (B5) + Laplace fast path (S6).
5. Statistical extensions: nugget (S2), Matérn (S1), diagnostics (S3), re-anchor (S4).
6. Spatio-temporal Kronecker-free path (B4).
7. Vignette, benchmarks, docs.
</content>
</invoke>
