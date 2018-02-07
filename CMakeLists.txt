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

include(GNUInstallDirs)
include(FeatureSummary)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(FDF REQUIRED libfdf>=0.1.1)
pkg_check_modules(PSOLVER REQUIRED psolver>=1.8)

include_directories(${FDF_INCLUDE_DIRS})
message(STATUS "Found FDF support: ${FDF_INCLUDE_DIRS};")
set(LIBS ${LIBS} ${FDF_LIBRARIES} ${PSOLVER_LIBRARIES})
link_directories(${FDF_LIBRARY_DIRS})

add_executable(${target_name}
  ${src_dir}/main.F90
  ${src_dir}/basis.F90
  ${src_dir}/density.F90
  ${src_dir}/hamiltonian.F90
  ${src_dir}/numeric.F90
  ${src_dir}/scf.F90
  ${src_dir}/smear.F90
  ${src_dir}/system.F90
  ${src_dir}/psolver.F90
  )
target_link_libraries(${target_name} ${FDF_LIBRARIES} ${PSOLVER_LIBRARIES})
target_include_directories(${target_name} PUBLIC ${FDF_INCLUDE_DIRS} ${PSOLVER_INCLUDE_DIRS})

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
