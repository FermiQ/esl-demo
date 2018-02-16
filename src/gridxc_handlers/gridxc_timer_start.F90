subroutine gridxc_timer_start(label)
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  character(len=*), intent(in) :: label

#if DEBUG_LEVEL > 1
  write(error_unit, '("[DEBUG] LibGridXC timer: START ",A)') label
#endif

end subroutine gridxc_timer_start
