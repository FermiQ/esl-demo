subroutine alloc_error_report(errmsg, errno)
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  character(len=*), intent(in) :: errmsg
  integer, intent(in)          :: errno

  write(error_unit, '("[DEBUG] LibGridXC memory error ", I12)') errno
  write(error_unit, '("          ",A)') trim(errmsg)

end subroutine alloc_error_report
