add_library(pw
  basis_pw.F90
  hamiltonian_pw.F90
  density_pw.F90
  utils_pw.F90
  )

target_include_directories(pw PRIVATE ${PROJECT_BINARY_DIR}/src/elsi_rci)
if(WITH_MPI)
  set_target_properties(pw PROPERTIES COMPILE_DEFINITIONS "WITH_MPI")
  target_link_libraries(pw elsi_rci common esldmpi)
else()
  target_link_libraries(pw elsi_rci common)
endif()
