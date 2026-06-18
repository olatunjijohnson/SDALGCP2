# SDALGCP2 — To-do

Drafted from a codebase scan on 2026-06-18. Items 1-4 completed 2026-06-19;
`R CMD check --as-cran` passes with 0 errors / 0 warnings / 4 NOTEs (all either
expected for a dev version or non-actionable — see below).

## Medium priority (API consistency)

- [ ] **Vignette build hygiene.** The precomputed `.png` figures the vignettes
  embed are copied into `inst/doc` and trip a `R CMD check` NOTE ("files should
  probably not be installed"). Decide: keep static figures (add a `.Rinstignore`
  / relocate sources), or compute figures live in the vignettes.

## Low priority / polish

- [ ] **Add a `LICENSE` file** (optional for GPL-2|GPL-3, but conventional).
- [ ] **Plan the path to `0.1.0`.** Current version is `0.0.0.9000`. Bumping to a
  proper `0.1.0` also clears the "Version contains large components" CRAN NOTE.

## Remaining R CMD check NOTEs (not blocking)

- *New submission / version `0.0.0.9000`* — expected for a dev version; clears at
  the `0.1.0` bump above.
- *`-mno-omit-leaf-frame-pointer` non-portable flag* — injected by the local R
  `Makeconf`, not our `Makevars`; not reproducible on CRAN.
- *PNGs installed under `inst/doc`* — the vignette-hygiene item above.

## Done

- [x] **Bundled example dataset `sdalgcp_data`** (64-region `sf`), documented,
  with a `data-raw/` generator and `LazyData: true`. (2026-06-19)
- [x] **Runnable `@examples`** added across the exported surface; all pass under
  `R CMD check --run-donttest`. (2026-06-19)
- [x] **Unified spatio-temporal prediction output**: `predict.SDALGCP2_ST` now
  returns a long `sf` with the same columns as the spatial predictor. (2026-06-19)
- [x] **R CMD check `--as-cran`** run and triaged: 0 errors, 0 warnings. (2026-06-19)
- [x] Prediction outputs renamed to `relative_risk`/`adjusted_rr`; `predict()`
  returns an `sf`. (commit 3dea05d)
- [x] Fixed continuous-prediction region assignment for points on shared
  polygon boundaries. (commit 3dea05d)
