module esl_info_m
  use iso_fortran_env, only : compiler_options,compiler_version
  use prec, only :ip 
  use yaml_output
#ifdef WITH_MPI
  use mpi
#endif
  implicit none
  private

  public :: about
contains

  subroutine about()

#ifdef WITH_MPI
    integer          :: ierr,mpi_ver,mpi_subver
    character(20)    :: aux
#ifndef OLDMPI    
    character(1000)  :: lib
    integer          :: ll
#endif    
#endif
    call yaml_mapping_open("About")
    call yaml_map("compiler version", trim(compiler_version()))
    call yaml_map("compiler options", trim(compiler_options()))
#ifdef WITH_MPI
    call MPI_Get_version(mpi_ver,mpi_subver,ierr)
    write(aux,'(i0,a1,i0)')mpi_ver,'.',mpi_subver
    call yaml_map("MPI Standard: ", trim(aux))
#ifndef OLDMPI    
    call MPI_Get_library_version(lib,ll,ierr)
    call yaml_map("MPI Implementation",trim(lib))
#endif    
#endif
    call yaml_mapping_close()
  end subroutine about

end module esl_info_m
