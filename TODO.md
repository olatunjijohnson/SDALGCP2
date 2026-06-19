# SDALGCP2 — To-do

The package is **CRAN-ready**: `R CMD check --as-cran` gives 0 errors, 0
warnings, and 2 NOTEs that are both expected (“New submission”; a
non-portable compiler flag injected by the local R `Makeconf`, not by
the package). Spelling is clean, all URLs resolve, every exported topic
has `\value`, and `cran-comments.md` documents the submission.

## Optional / nice-to-have (not blocking submission)

**Convert the remaining vignettes to live evaluation** (raster,
misaligned, spatio-temporal, scale, confounding). They currently stay
`eval = FALSE` with precomputed figures (kept out of the install via
`.Rinstignore`), which is CRAN-clean and keeps build times low. Making
them live would be more reproducible but slower to build.

**A `liver` case-study vignette** — the real dataset end to end.

Add the maintainer’s ORCID to `Authors@R` if desired.

## Done

**Made the package fully independent of SDALGCP**: the `liver` dataset
builds from a frozen `data-raw/liver_source.rds` (no SDALGCP needed);
removed all comparative framing from
DESCRIPTION/roxygen/README/NEWS/vignette/comments/design docs; removed
the `*_vs_SDALGCP` scripts. Method paper citations retained.
(2026-06-19)

**CRAN prep**: `.Rinstignore` for the `t2`-`t6` vignette PNGs (clears
the inst/doc NOTE); `Language: en-GB` + `inst/WORDLIST` (spelling
clean); `cran-comments.md`; final `--as-cran` run clean. (2026-06-19)

**Version 0.1.0** + build-ignored `LICENSE.md` (GPL-2 \| GPL-3).
(2026-06-19)

**Real dataset `liver`** (PBC by LSOA, Johnson et al. 2019) with
attribution; `data-raw/liver.R`. (2026-06-19)

**Intro vignette runs live on `sdalgcp_data`**; `t1_*` figures retired.
(2026-06-19)

**Simulated dataset `sdalgcp_data`** + `LazyData: true`. (2026-06-19)

**Runnable `@examples`** across the exported surface. (2026-06-19)

**Unified spatio-temporal prediction output** to a long `sf`.
(2026-06-19)

Prediction outputs renamed to `relative_risk`/`adjusted_rr`;
[`predict()`](https://rdrr.io/r/stats/predict.html) returns an `sf`;
continuous-prediction boundary-point fix. (commit 3dea05d)
