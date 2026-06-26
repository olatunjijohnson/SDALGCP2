library(testthat)
library(SDALGCP2)

# Cap parallelism for CRAN's check machines (policy: use at most 2 cores).
# A runtime omp_set_num_threads() call is used rather than the OMP_NUM_THREADS
# environment variable: on some platforms the OpenMP runtime is already
# initialised by the time the tests run, so the variable is read too late and
# the parallel C++ routines otherwise use every core (CPU time >> elapsed time).
SDALGCP2:::set_omp_num_threads(2L)

test_check("SDALGCP2")
