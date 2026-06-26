# Cap parallelism for CRAN's check machines (policy: use at most 2 cores).
# Set before SDALGCP2's shared library is loaded so libgomp reads the limit
# at initialisation; otherwise the OpenMP regions in the C++ code default to
# all available cores and CPU time greatly exceeds elapsed time.
Sys.setenv(OMP_NUM_THREADS = "2")
Sys.setenv(OMP_THREAD_LIMIT = "2")

library(testthat)
library(SDALGCP2)

test_check("SDALGCP2")
