## Submission

This is the first CRAN submission of SDALGCP2 (version 0.1.0).

SDALGCP2 fits a spatially discrete approximation to a log-Gaussian Cox process for
aggregated disease count data, with the performance-critical steps implemented in
C++ (RcppArmadillo).

## Test environments

* local: Ubuntu 24.04, R 4.5.2
* GitHub Actions (r-lib/actions, `--as-cran`): Ubuntu, Windows, macOS, R release

## R CMD check results

0 errors | 0 warnings | 2 NOTEs

* "New submission." This is the first submission of the package.

* "Compilation used the following non-portable flag(s): `-mno-omit-leaf-frame-pointer`."
  This flag is injected by the local R installation's `Makeconf`, not by the
  package's `src/Makevars` (which sets only the OpenMP flags). It does not appear
  on the standard CRAN build machines.

The compiled shared library is reported as ~3.9 Mb (INFO, not a NOTE) because the
C++ kernels use RcppArmadillo templates; there is no R-level way to reduce it.

## Reverse dependencies

There are no reverse dependencies (new package).

## Data

The bundled `liver` dataset is the disease-count study data of Johnson, Diggle
and Giorgi (2019) <doi:10.1002/sim.8339> (authored by the package maintainer),
aggregated to LSOA level; the source is attributed in the dataset's
documentation. The other bundled dataset, `sdalgcp_data`, is simulated.
