add_library(pw
  basis_pw.F90
  density_pw.F90
  )

target_link_libraries(pw common)
