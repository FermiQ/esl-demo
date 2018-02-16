subroutine alloc_memory_event(bytes, label)
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  integer, intent(in)          :: bytes
  character(len=*), intent(in) :: label

#if DEBUG_LEVEL > 1
  write(error_unit, '("[DEBUG] LibGridXC allocated ",F8.3, "Mbytes (",A,")")') &
&   bytes, label
#endif

end subroutine alloc_memory_event
