option(WITH_MPI "build a MPI version" ON)
option(WITH_FLOOK "Add flook library support to interface Lua scripts" ON)
option(WITH_OPENMP "build an OpenMP version" OFF)
option(BUILD_TESTING "Build with Testing support" OFF)
option(WITH_COVERAGE "Build with instrumentation for code coverage" OFF)
option(WITH_DOC "Doxygen Documentation" OFF)
option(BUILD_SHARED_LIBS "Build with shared libraries" OFF)

set(MPI_NPROCS 8 CACHE STRING "number of MPI processes to be used for code coverage and tests")
