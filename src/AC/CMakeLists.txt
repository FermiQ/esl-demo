add_library(ac
  basis_ac.F90
  create_sparse_pattern.F90
  density_matrix.F90  
  mulliken.F90
  overlap_matrix.F90
  sparse_pattern.F90
  sparse_matrix.F90
  )

target_link_libraries(ac common)
