# SDALGCP2 — To-do

Most release-readiness work is done (see below). `R CMD check --as-cran` is clean
apart from expected/non-actionable NOTEs.

## Remaining

- [ ] **Convert the remaining vignettes to live evaluation** (raster, misaligned,
  spatio-temporal, scale, confounding). The intro vignette is now live and
  figure-free; the others still embed precomputed `t2_*`–`t6_*` PNGs, which is the
  sole remaining `R CMD check` NOTE ("files should probably not be installed").
  Converting them (or adding a real case study on `liver`) clears it for good.
- [ ] **A `liver` case-study vignette** — use the real bundled dataset end to end
  (deprivation covariates, exceedance hotspots) as a realistic companion to the
  simulated intro.

## Remaining R CMD check NOTEs (not blocking)

- *`-mno-omit-leaf-frame-pointer` non-portable flag* — injected by the local R
  `Makeconf`, not our `Makevars`; not reproducible on CRAN.
- *PNGs installed under `inst/doc`* — the vignette-conversion item above (now only
  `t2`–`t6`, since the intro no longer ships figures).

## Done

- [x] **Version 0.1.0** and a build-ignored `LICENSE.md` (GPL-2 | GPL-3 notice);
  the version bump clears the "large components" CRAN NOTE. (2026-06-19)
- [x] **Real dataset `liver`** — primary biliary cirrhosis counts by LSOA, derived
  with attribution from `SDALGCP::PBCshp_sf` (Johnson et al. 2019); see
  `data-raw/liver.R`. (2026-06-19)
- [x] **Intro vignette runs live on `sdalgcp_data`**; `t1_*` figures retired. (2026-06-19)
- [x] **Bundled simulated dataset `sdalgcp_data`** + `data-raw/` generator +
  `LazyData: true`. (2026-06-19)
- [x] **Runnable `@examples`** across the exported surface; pass under
  `R CMD check --run-donttest`. (2026-06-19)
- [x] **Unified spatio-temporal prediction output** to a long `sf`. (2026-06-19)
- [x] Prediction outputs renamed to `relative_risk`/`adjusted_rr`; `predict()`
  returns an `sf`; continuous-prediction boundary-point fix. (commit 3dea05d)
