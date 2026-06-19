## Submission

This is the first CRAN submission of SDALGCP2 (version 0.1.0).

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
