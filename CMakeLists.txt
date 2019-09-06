cmake_minimum_required(VERSION 3.02)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0048 NEW)
endif(COMMAND cmake_policy)

project(esl-demo
  VERSION 0.0.1
  LANGUAGES Fortran CXX
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
if (WITH_COVERAGE)
  include(CodeCoverage)
  set(BUILD_TESTING ON)
endif() 

if(WITH_MPI)
  ### if we want MPI check for it and set the internal kitchen
  find_package(MPI REQUIRED)
  add_definitions(${MPI_Fortran_COMPILE_FLAGS})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_Fortran_LIBRARIES})
  message(STATUS "Build with MPI support!")
  add_definitions("-DWITH_MPI")
else()
  message(STATUS "Build without MPI support!")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(FDF REQUIRED libfdf>=0.1.1)
pkg_check_modules(PSOLVER REQUIRED psolver>=1.8)
pkg_check_modules(PSPIO REQUIRED libpspio>=0.2.2)
pkg_check_modules(FFTW REQUIRED fftw3>=3.3.0)
pkg_check_modules(GSL REQUIRED gsl>=1.16)
pkg_check_modules(GRIDXC REQUIRED libgridxc>=0.8.0.3)

message(STATUS "Found FDF support!")
include_directories(${FDF_INCLUDE_DIRS})
link_directories(${FDF_LIBRARY_DIRS})

message(STATUS "Found PSOLVER support!")
include_directories(${PSOLVER_INCLUDE_DIRS})
link_directories(${PSOLVER_LIBRARY_DIRS})

message(STATUS "Found PSPIO support!")
include_directories(${PSPIO_INCLUDE_DIRS})
link_directories(${PSPIO_LIBRARY_DIRS})

message(STATUS "Found GSL support!")
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

message(STATUS "Found FFTW support!")
include_directories(${FFTW_PREFIX}/include)
link_directories(${FFTW_LIBRARY_DIRS})

message(STATUS "Found GRIDXC support!")
include_directories(${GRIDXC_INCLUDEDIR})
include_directories(${GRIDXC_INCLUDE_DIRS})
link_directories(${GRIDXC_LIBRARY_DIRS})

if(WITH_ELSI)
  pkg_check_modules(ELSI REQUIRED elsi)
else()
  pkg_check_modules(ELSI elsi)
endif()
if(ELSI_FOUND)
  if(NOT DEFINED WITH_ELSI)
    set(WITH_ELSI true)
    message(STATUS "ELSI will be enabled by default!")
  endif()
endif()
if(WITH_ELSI)
  include_directories(${ELSI_INCLUDEDIR})
  include_directories(${ELSI_INCLUDE_DIRS})
  link_directories(${ELSI_LIBRARY_DIRS})
  message(STATUS "Build with ELSI library support!")
  add_definitions("-DWITH_ELSI")
else()
  unset(ELSI_INCLUDEDIR)
  unset(ELSI_INCLUDE_DIRS)
  unset(ELSI_LIBRARY_DIRS)
  unset(ELSI_LIBRARIES)
  message(STATUS "Build without ELSI library, no advanced Eigensolver support (PW won't work)!")
endif()

if(WITH_FLOOK)
  pkg_check_modules(FLOOK REQUIRED flook)
  include_directories(${FLOOK_INCLUDE_DIRS})
  include_directories(${FLOOK_INCLUDEDIR})
  link_directories(${FLOOK_LIBRARY_DIRS})
  message(STATUS "Build with flook library and Lua-interface support!")
  add_definitions("-DWITH_FLOOK")

  # Also add requirement of fdict
  pkg_check_modules(FDICT REQUIRED fdict)
  include_directories(${FDICT_INCLUDEDIR})
  include_directories(${FDICT_INCLUDE_DIRS})
  link_directories(${FDICT_LIBRARY_DIRS})

else()
  message(STATUS "Build without flook library (no Lua-interface support)!")
endif()


include_directories(${PROJECT_BINARY_DIR}/src/common)
include_directories(${PROJECT_BINARY_DIR}/src/AC)
include_directories(${PROJECT_BINARY_DIR}/src/PW)

add_subdirectory(${src_dir}/PW)
add_subdirectory(${src_dir}/AC)
add_subdirectory(${src_dir}/common)
add_subdirectory(${src_dir}/elsi_rci)
add_subdirectory(${src_dir}/gridxc_handlers)


add_executable(${target_name}
  ${src_dir}/main/basis.F90
  ${src_dir}/main/density.F90
  ${src_dir}/main/hamiltonian.F90
  ${src_dir}/main/main.F90
  ${src_dir}/main/next_step.F90
  ${src_dir}/main/scf.F90
  ${src_dir}/main/system.F90
  )
add_dependencies(${target_name} gridxc_handlers)
target_link_libraries(${target_name} ac pw common ${FDF_LIBRARIES} ${PSOLVER_LIBRARIES} ${PSPIO_LIBRARIES} ${GRIDXC_LIBRARIES} gridxc_handlers ${FFTW_LIBRARIES} ${GSL_LIBRARIES})
if(WITH_MPI)
  set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "WITH_MPI")
  target_link_libraries(${target_name} esldmpi ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()
if(WITH_ELSI)
  set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "WITH_ELSI")
  target_link_libraries(${target_name} ${ELSI_LIBRARIES})
endif()
if(WITH_FLOOK)
  set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "WITH_FLOOK")
  target_link_libraries(${target_name}  ${FLOOK_LIBRARIES} flookall ${FDICT_LIBRARIES})
endif()

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

if(BUILD_TESTING)
  include(CTest)
  message(STATUS "Build with testing support support: make test")
  find_package(PythonInterp 3.0 REQUIRED)
  macro (do_test arg)
    if (EXISTS ${CMAKE_SOURCE_DIR}/tests/${arg}/${arg}.inp)
      message(STATUS "Adding test ${arg}/${arg}.inp")
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/${arg})
      file(COPY ${CMAKE_SOURCE_DIR}/tests/${arg}/
        DESTINATION ${CMAKE_BINARY_DIR}/Testing/${arg}/)
      add_test(NAME ${arg}
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/utils/checktest.py -c ${CMAKE_SOURCE_DIR}/tests/esl-testsuite.yml 
          -r ${CMAKE_SOURCE_DIR}/utils/Refs -s ${arg} -e ${CMAKE_BINARY_DIR}/bin/${target_name} 
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/${arg})
        set_tests_properties (${arg}
        PROPERTIES PASS_REGULAR_EXPRESSION "Status: PASSED")
    endif()
  endmacro (do_test)
  file(GLOB Ts LIST_DIRECTORIES true ${CMAKE_SOURCE_DIR}/tests/*)
  foreach(T ${Ts})
    string(REGEX MATCH "([a-zA-Z0-9-_]+)$" test ${T})
    do_test(${test})
  endforeach()
else()
  message(STATUS "Build without testing target")
endif()


if (WITH_COVERAGE)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/analysis)
  setup_coverage(coverage analysis)
endif()

#installation
install(TARGETS ${target_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
