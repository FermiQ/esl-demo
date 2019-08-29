set(common_SOURCES
  basis_base.F90
  constants.F90
  elsi_wrapper.F90
  energy.F90
  esl_dict.F90
  flook_global.F90
  flook_if.F90
  force.F90
  geometry.F90
  grid.F90
  info.F90
  ion_interaction.F90
  message.F90
  mixing.F90
  numeric.F90
  potential.F90
  prec.F90
  psolver.F90
  smear.F90
  species.F90
  states.F90
  xc.F90  
  )

if(NOT WITH_ELSI)
  list(APPEND common_SOURCES elsi_fake.F90)
endif()

add_library(common ${common_SOURCES})

if(WITH_MPI)
  add_library(esldmpi
    mpi_dist_block_cyclic.F90 
    mpi_dist_cyclic.F90
    mpi_dist.F90
    )
  set_target_properties(esldmpi PROPERTIES
    COMPILE_DEFINITIONS "WITH_MPI"
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    )
  set_target_properties(common PROPERTIES
    COMPILE_DEFINITIONS "WITH_MPI"
    )
endif()
if (WITH_FLOOK)
  set_target_properties(common PROPERTIES COMPILE_DEFINITIONS "WITH_FLOOK")
endif()
