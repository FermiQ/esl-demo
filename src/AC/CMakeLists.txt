add_library(ac
  basis_ac.F90
  basis_ac_abc.F90
  calc_density_matrix.F90
  create_sparse_pattern.F90
  density_matrix.F90  
  hamiltonian_ac.F90  
  mulliken_ac.F90
  density_ac.F90
  overlap_matrix.F90
  sparse_pattern.F90
  sparse_matrix.F90
  )

target_link_libraries(ac common)
if(WITH_MPI)
  set_target_properties(ac PROPERTIES COMPILE_DEFINITIONS "WITH_MPI" )
endif()
