module esl_flook_global_m

#ifdef WITH_FLOOK
  use flook, only : luaState
#endif

  implicit none

#ifdef WITH_FLOOK
  ! LUA-handle, we really do need it up here to be able to use it
  type(luaState) :: LUA
#else
  integer, parameter :: LUA_NEVER_USED = 0!!
#endif

end module esl_flook_global_m


