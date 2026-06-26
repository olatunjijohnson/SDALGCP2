## Submission

This is a resubmission of the first CRAN submission of SDALGCP2 (version 0.1.0).
A previous upload was auto-rejected at the incoming pre-test. This version
addresses every issue raised, including the Debian tests NOTE that persisted on
the most recent pre-test:

* The two README links to `math/*.pdf` (flagged as invalid file URIs) now use
  absolute GitHub URLs; those PDFs live in the repository's `math/` directory,
  which is excluded from the package build via `.Rbuildignore`.

* "Running R code in 'testthat.R' had CPU time 5.x times elapsed time" (Debian):
  the parallel C++ routines use OpenMP, and on the pre-test machine the runtime
  was already initialised before the tests ran, so setting `OMP_NUM_THREADS`
  from R had no effect. The test runner now caps threads with a *runtime*
  `omp_set_num_threads(2)` call (an internal helper), which takes effect
  regardless of when the OpenMP runtime was initialised. End users are
  unaffected: the exported functions still default to the OpenMP default and
  expose an `nthreads` argument.

* The DESCRIPTION wording was adjusted so the spell checker no longer flags
  words: the method names/acronyms `'Matern'`, `'MALA'` and `'MCML'` are quoted,
  the citation is written in full ("Johnson, Diggle and Giorgi (2019)"), and
  "spatio-temporal" is written as "space-time".

SDALGCP2 fits a spatially discrete approximation to a log-Gaussian Cox process for
aggregated disease count data, with the performance-critical steps implemented in
C++ (RcppArmadillo).

## Test environments

* local: Ubuntu 24.04, R 4.5.2 (`R CMD check --as-cran`)
* win-builder: R-release (R 4.6.0) and R-devel
* macOS builder (mac.r-project.org): R-release
* GitHub Actions (r-lib/actions): macOS, Windows, and Ubuntu (R-devel, release, oldrel)

## R CMD check results

0 errors | 0 warnings | 1 note

* "New submission." This is the first submission of the package.

The installed size is ~5.9 Mb because the compiled shared library uses
RcppArmadillo templates (reported as INFO, not a NOTE).

(On the maintainer's local machine a second NOTE about a non-portable compiler
flag, `-mno-omit-leaf-frame-pointer`, also appears; this flag is injected by the
local R installation's `Makeconf`, not by the package's `src/Makevars`, and does
not occur on win-builder or the other test environments.)

## Reverse dependencies

There are no reverse dependencies (new package).

## Data

The bundled `liver` dataset is the disease-count study data of Johnson, Diggle
and Giorgi (2019) <doi:10.1002/sim.8339> (authored by the package maintainer),
aggregated to LSOA level; the source is attributed in the dataset's
documentation. The other bundled dataset, `sdalgcp_data`, is simulated.
