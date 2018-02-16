subroutine die(errmsg)
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  character(len=*), intent(in) :: errmsg

  write(error_unit, '(A)') "[DEBUG] LibGridXC error:"
  write(error_unit, '("          ",A)') trim(errmsg)

  stop

end subroutine die
