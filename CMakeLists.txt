cmake_minimum_required(VERSION 3.02)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0048 NEW)
endif(COMMAND cmake_policy)

project(esl-demo
  VERSION 0.0.1
  LANGUAGES Fortran
  )

set(AUTHOR "X;Y;Z")
set(AUTHOR_DETAILS "")
set(DESCRIPTION "ESL Demonstator")

set(src_dir ${CMAKE_SOURCE_DIR}/src)
set(target_name esl-demo.X)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include(FeatureSummary)
include(ESLDemoBuildOptions)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

if(WITH_MPI)
  ### if we want MPI check for it and set the internal kitchen
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "Build with MPI support!")
else()
  message(STATUS "Build without MPI support!")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(FDF REQUIRED libfdf>=0.1.1)
pkg_check_modules(PSOLVER REQUIRED psolver>=1.8)

include_directories(${FDF_INCLUDE_DIRS})
message(STATUS "Found FDF support: ${FDF_INCLUDE_DIRS};")
link_directories(${FDF_LIBRARY_DIRS})


if(WITH_MPI)
  add_library(esldmpi
    ${src_dir}/mpi_dist_block_cyclic.F90 
    ${src_dir}/mpi_dist_cyclic.F90
    ${src_dir}/mpi_dist.F90
    )
  set_target_properties(esldmpi PROPERTIES
    COMPILE_DEFINITIONS "WITH_MPI"
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    )
endif()


add_executable(${target_name}
  ${src_dir}/main.F90
  ${src_dir}/basis.F90
  ${src_dir}/density.F90
  ${src_dir}/energy.F90
  ${src_dir}/grid.F90
  ${src_dir}/hamiltonian.F90
  ${src_dir}/init_sparse_pattern.F90
  ${src_dir}/message.F90
  ${src_dir}/mixing.F90
  ${src_dir}/new_step.F90
  ${src_dir}/numeric.F90
  ${src_dir}/psolver.F90
  ${src_dir}/potential.F90
  ${src_dir}/scf.F90
  ${src_dir}/smear.F90
  ${src_dir}/sparse_pattern.F90
  ${src_dir}/sparse_matrix.F90
  ${src_dir}/states.F90
  ${src_dir}/system.F90
  ${src_dir}/elsi_wrapper.F90
  ${src_dir}/elsi_fake.F90
  ${src_dir}/xc.F90
  )
target_link_libraries(${target_name} ${FDF_LIBRARIES} ${PSOLVER_LIBRARIES} )
if(WITH_MPI)
  target_link_libraries(${target_name} esldmpi ${MPI_Fortran_LIBRARIES})
endif()
target_include_directories(${target_name} PUBLIC ${FDF_INCLUDE_DIRS} ${PSOLVER_INCLUDE_DIRS} )

if(WITH_DOC)
  message(STATUS "Build Doxygen API Documentation")
  find_package(Doxygen REQUIRED)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/Doxyfile.cmake Doxyfile)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile)
  install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} doc)")
  install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/developers/)
endif()

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
