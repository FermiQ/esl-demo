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
pkg_check_modules(PSPIO REQUIRED libpspio>=0.2.2)

include_directories(${FDF_INCLUDE_DIRS})
message(STATUS "Found FDF support: ${FDF_INCLUDE_DIRS};")
link_directories(${FDF_LIBRARY_DIRS})

include_directories(${PSPIO_INCLUDE_DIRS})
message(STATUS "Found PSPIO support: ${PSPIO_INCLUDE_DIRS};")
link_directories(${PSPIO_LIBRARY_DIRS})


if(WITH_FLOOK)
  pkg_check_modules(FLOOK REQUIRED flook)
  include_directories(${FLOOK_Fortran_INCLUDE_PATH})
  message(STATUS "Build with flook library and Lua-interface support!")
  add_definitions("-DWITH_FLOOK")
else()
  message(STATUS "Build without flook library (no Lua-interface support)!")
endif()


include_directories(${PROJECT_BINARY_DIR}/src/common)
include_directories(${PROJECT_BINARY_DIR}/src/AC)
include_directories(${PROJECT_BINARY_DIR}/src/PW)

add_subdirectory(${src_dir}/PW)
add_subdirectory(${src_dir}/AC)
add_subdirectory(${src_dir}/common)


add_executable(${target_name}
  ${src_dir}/main/basis.F90
  ${src_dir}/main/density.F90
  ${src_dir}/main/hamiltonian.F90
  ${src_dir}/main/main.F90
  ${src_dir}/main/next_step.F90
  ${src_dir}/main/scf.F90
  ${src_dir}/main/system.F90
  )
target_link_libraries(${target_name} ac pw common ${FDF_LIBRARIES} ${PSOLVER_LIBRARIES} ${PSPIO_LIBRARIES})
if(WITH_MPI)
  target_link_libraries(${target_name} esldmpi ${MPI_Fortran_LIBRARIES})
endif()
target_include_directories(${target_name} PUBLIC ${FDF_INCLUDE_DIRS} ${PSOLVER_INCLUDE_DIRS} ${PSPIO_INCLUDE_DIRS})

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

