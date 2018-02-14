add_library(pw
  basis_pw.F90
  hamiltonian_pw.F90
  density_pw.F90
  utils_pw.F90
  )

if(WITH_MPI)
  target_link_libraries(pw elsi_rci common esldmpi)
else()
  target_link_libraries(pw elsi_rci common)
endif()
