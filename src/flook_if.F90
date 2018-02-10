!< Module to interface flook library to ESL-demo.
!<
!< The flook library enables running Lua scripts from within
!< the fortran code.
!< This code is heavily interacting with esl_dicts_m module
!< and the dictionaries within.
module esl_flook_if_m

#ifdef WITH_LUA
  use flook
#endif

  implicit none

  private

  ! Signals to LUA
  ! Right after reading initial options 
  integer, parameter, public :: LUA_INITIALIZE = 1
  ! Right before SCF step starts, but at each step
  integer, parameter, public :: LUA_INIT_STEP = 2
  ! at the start of each SCF step
  integer, parameter, public :: LUA_SCF_LOOP = 3
  ! after each SCF has finished
  integer, parameter, public :: LUA_FORCES = 4
  ! when moving the atoms, right after the FORCES step
  integer, parameter, public :: LUA_NEXT_STEP = 5

#ifdef WITH_LUA

  public :: flook_if_init, flook_if_call, flook_if_close

  ! Internal parameters
  logical, save :: flook_if_run = .false.
  character(len=512), save, public :: flook_if_file = ' '
  ! Debugging flag for both parallel and serial debugging
  logical, save, public :: flook_if_debug = .false.

contains

  subroutine flook_if_init(LUA)

    use fdf, only : fdf_get

    type(luaState), intent(inout) :: LUA

    character(len=30) :: fortran_msg

    character(*), parameter :: fortran_static_lua = '&
esl = { &
  Node = 1, &
  INITIALIZE = 1, &
  INIT_STEP = 2, &
  SCF_LOOP = 3, &
  FORCES = 4, &
  NEXT_STEP = 5, &
  ANALYSIS = 6, &
  state = 0, &
  IOprint = function(self, ...) &
    if self.IONode then &
       print(...) &
     end &
  end, &
  print = function(self, ...) &
    print(...) &
  end, &
} &
IOprint = function(...) &
  esl:IOprint(...) &
end &
esl_comm = function(...) end'

    ! For error-handling with lua
    integer :: err
    character(len=2048) :: err_msg

    ! First retrieve lua file
    flook_if_file = fdf_get('LUA.Script',' ')
    ! Immediately return if the file is not specified...
    if ( len_trim(flook_if_file) == 0 ) return

    ! Default debugging only on the io-node.
    flook_if_debug = fdf_get('LUA.Debug',.false.)
    if ( fdf_get('LUA.Debug.MPI',.false.) ) then
      ! Only if requesting parallel debug should all processors
      ! use the debugging.
      flook_if_debug = .true.
    end if

    ! Initialize the Lua state
    call lua_init(LUA)

    ! Create LUA table for data container
    call lua_run(LUA, code = fortran_static_lua )

    call lua_register(LUA,'esl_receive', flook_if_receive)
    call lua_register(LUA,'esl_send', flook_if_send)
    ! Make local esl.receive and esl.send
    call lua_run(LUA, code = 'esl.receive = esl_receive' )
    call lua_run(LUA, code = 'esl.send = esl_send' )

    ! Only used for printing information about
    ! what can be retrieved
    ! This function will return different things
    ! dependen on where in the routine it is called
    call lua_register(LUA,'_internal_print_allowed', flook_if_print_objects)
    call lua_run(LUA, code = 'esl.print_allowed = _internal_print_allowed' )

    ! TODO when adding parallelism
    write(fortran_msg,'(a,i0)') 'esl.Node = ', 1
    call lua_run(LUA, code = fortran_msg )
    write(fortran_msg,'(a,i0)') 'esl.Nodes = ',1
    call lua_run(LUA, code = fortran_msg )
    call lua_run(LUA, code = 'esl.IONode = true' )

    ! Run the requested lua-script
    err_msg = " "
    call lua_run(LUA, flook_if_file, error = err, message=err_msg)
    if ( err /= 0 ) then
      write(*,'(a)') trim(err_msg)
      call message_error('LUA initialization failed, please check your Lua script!!!')
    end if

  end subroutine flook_if_init

  subroutine flook_if_call(LUA, state)
    type(luaState), intent(inout) :: LUA
    integer, intent(in) :: state
    character(len=30) :: tmp

    ! For error-handling with lua
    integer :: err
    character(len=2048) :: err_msg

    ! Return immediately if we should not run
    if ( .not. flook_if_run ) return

    ! Transfer the state to the lua interpreter such
    ! that decisions can be made as to which steps
    ! to take
    write(tmp,'(a,i0)') 'esl.state = ',state
    call lua_run(LUA, code = tmp )

    if ( flook_if_debug ) then
      write(*,'(a,i0)') 'esl-lua: calling esl_comm() @ ',state
    end if

    ! Call communicator
    call lua_run(LUA, code = 'esl_comm()', error = err, message=err_msg )
    if ( err /= 0 ) then
      write(*,'(a)') trim(err_msg)
      call message_error('LUA could not run esl_comm() without an error, please &
          &check your Lua script')
    end if

  end subroutine flook_if_call

  subroutine flook_if_close(LUA)
    type(luaState), intent(inout) :: LUA
    
    ! Return immediately if we should not run
    if ( .not. flook_if_run ) return
    call lua_close(LUA)

  end subroutine flook_if_close


  ! ! ! ! ! ! ! 
  ! The remaining functions/routines are private
  ! methods.
  ! ! ! ! ! ! !


  function flook_if_receive(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use esl_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl

    type(dict) :: keys

    if ( flook_if_debug ) then
      write(*,'(a,i0)') '  lua: esl_receive, Node = ',Node + 1
    end if

    call lua_init(LUA,state)

    ! Retrieve information
    call flook_if_get_tbl_to_dict(lua,keys)

    ! open global esl table
    tbl = lua_table(LUA,'esl')

    ! Expose the dictionary
    call flook_if_put_dict(tbl,options,keys)
    call flook_if_put_dict(tbl,variables,keys)

    call lua_close_tree(tbl)

    ! Clean-up
    call delete(keys)

    ! this function returns nothing
    nret = 0

  end function flook_if_receive

  function flook_if_send(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use esl_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl
    type(dict) :: keys

    if ( flook_if_debug ) then
      write(*,'(a,i0)') '  lua: esl_send, Node = ',Node + 1
    end if

    call lua_init(LUA,state)

    ! Retrieve information
    call flook_if_get_tbl_to_dict(lua,keys)

    ! open global esl table
    tbl = lua_table(LUA,'esl')

    ! Expose the dictionary
    call flook_if_get_dict(tbl,options,keys)
    call flook_if_get_dict(tbl,variables,keys)

    call lua_close_tree(tbl)

    ! Clean-up
    call delete(keys)

    ! this function returns nothing
    nret = 0

  end function flook_if_send

  subroutine flook_if_get_tbl_to_dict(lua,keys)

    use variable
    use dictionary

    type(luaState), intent(inout) :: lua
    type(dict), intent(inout) :: keys

    type(luaTbl) :: tbl
    character(len=255) :: name
    integer :: i, N

    ! Clean the dictionary
    call delete(keys)

    ! Retrieve the table @ the top
    tbl = lua_table(LUA)

    ! Traverse all elements in the table
    N = len(tbl)
    do i = 1 , N
      call lua_get(tbl,i,name)
      keys = keys // (trim(name).kv.1)
    end do
    ! Loop through all keys in it.
    !    print *,'Number of elements passed: ',len(tbl),len(keys)
    !    call print(keys)

    call lua_close(tbl)

  end subroutine flook_if_get_tbl_to_dict

  subroutine flook_if_put_dict(tbl,dic,keys)

    use variable
    use dictionary

    type(luaTbl), intent(inout) :: tbl
    type(dict), intent(inout) :: dic
    type(dict), intent(inout), optional :: keys

    ! Sadly we need a pointer for all variables that we might
    ! expect to handle.
    character(len=DICT_KEY_LENGTH) :: key
    type(dict) :: pd ! pointer to the dictionary
    type(var) :: v

    if ( present(keys) ) then

      ! Loop over all entries in the keys dictionary
      pd = .first. keys
      do while ( .not. (.empty. pd) )
        key = trim(.key. pd)
        if ( key .in. dic ) then
          call associate(v,dic,key)
          call put_var(key)
          call nullify(v) ! do not delete, simply nullify
        end if
        pd = .next. pd
      end do

    else

      ! Loop over all entries
      pd = .first. dic
      do while ( .not. (.empty. pd) )
        key = .key. pd
        call associate(v,dic,trim(key))
        call put_var(key)
        call nullify(v)
        pd = .next. pd
      end do

    end if

  contains

    subroutine put_var(key)
      use prec, only: sp, dp
      character(len=*), intent(in) :: key
      character(len=1), pointer :: a1(:)
      logical, pointer :: b0, b1(:), b2(:,:)
      integer, pointer :: i0, i1(:), i2(:,:)
      real(sp), pointer :: s0, s1(:), s2(:,:)
      real(dp), pointer :: d0, d1(:), d2(:,:)
      character(len=4) :: t
      character(len=255) :: lkey, rkey
      integer :: lvls
      lvls = 0
      t = which(v)
      if ( flook_if_debug ) then
        write(*,'(4a)') '    esl2lua; dtype = ',t,', var = ',trim(key)
      end if
      !      print *,'Attempt storing: ',trim(key), ' type= ',t
      select case ( t )
      case ( 'a1', 'b0' , 'i0', 's0', 'd0' )
        ! We need to handle the variables secondly
        call key_split(key,lkey,rkey)
        if ( len_trim(rkey) == 0 ) then
          lvls = 0
          rkey = lkey
        else
          call lua_open(tbl,lkey,lvls = lvls)
        end if
      end select
      select case ( t )
      case ( 'a1' )
        call associate(a1,v)
        lkey = cunpack(a1)
        call lua_set(tbl,rkey,lkey(1:len_trim(lkey)))
      case ( 'b0' )
        call associate(b0,v)
        call lua_set(tbl,rkey,b0)
      case ( 'b1' )
        call associate(b1,v)
        call lua_set(tbl,key,b1)
      case ( 'b2' ) 
        call associate(b2,v)
        call lua_set(tbl,key,b2)
      case ( 'i0' ) 
        call associate(i0,v)
        call lua_set(tbl,rkey,i0)
      case ( 'i1' ) 
        call associate(i1,v)
        call lua_set(tbl,key,i1)
      case ( 'i2' ) 
        call associate(i2,v)
        call lua_set(tbl,key,i2)
      case ( 's0' ) 
        call associate(s0,v)
        call lua_set(tbl,rkey,s0)
      case ( 's1' ) 
        call associate(s1,v)
        call lua_set(tbl,key,s1)
      case ( 's2' ) 
        call associate(s2,v)
        call lua_set(tbl,key,s2)
      case ( 'd0' ) 
        call associate(d0,v)
        call lua_set(tbl,rkey,d0)
        !         print *,'setting: '//trim(key)//' to ',d0
      case ( 'd1' ) 
        call associate(d1,v)
        call lua_set(tbl,key,d1)
      case ( 'd2' ) 
        call associate(d2,v)
        call lua_set(tbl,key,d2)
        !         print *,'setting: '//trim(key)//' to ',d2
      end select
      !      print *,'Done storing: ',trim(key), ' type= ',t
      call lua_close(tbl,lvls = lvls)
    end subroutine put_var

  end subroutine flook_if_put_dict

  subroutine flook_if_get_dict(tbl,dic,keys)

    use variable
    use dictionary

    type(luaTbl), intent(inout) :: tbl
    type(dict), intent(inout) :: dic
    type(dict), intent(in), optional :: keys

    ! Sadly we need a pointer for all variables that we might
    ! expect to handle.
    character(len=DICT_KEY_LENGTH) :: key
    type(dict) :: pd ! pointer to the dictionary
    type(var) :: v

    if ( present(keys) ) then

      ! Loop over all entries in the keys dictionary
      pd = .first. keys
      do while ( .not. (.empty. pd) )
        key = .key. pd
        if ( key .in. dic ) then
          call associate(v,dic,trim(key))
          call get_var(key)
          call nullify(v)
        end if
        pd = .next. pd
      end do

    else

      ! Loop over all entries
      pd = .first. dic
      do while ( .not. (.empty. pd) )
        key = .key. pd
        call associate(v,dic,trim(key))
        call get_var(key)
        call nullify(v)
        pd = .next. pd
      end do

    end if

  contains

    subroutine get_var(key)
      use prec, only: sp, dp
      character(len=*), intent(in) :: key
      character(len=256) :: V0
      character(len=1), pointer :: a1(:)
      logical, pointer :: b0, b1(:), b2(:,:)
      integer, pointer :: i0, i1(:), i2(:,:)
      real(sp), pointer :: s0, s1(:), s2(:,:)
      real(dp), pointer :: d0, d1(:), d2(:,:)
      character(len=4) :: t
      character(len=255) :: lkey, rkey
      integer :: na1
      integer :: lvls
      lvls = 0
      t = which(v)
      if ( flook_if_debug ) then
        write(*,'(4a)') '    lua2esl; dtype = ',t,', var = ',trim(key)
      end if
      select case ( t )
      case ( 'a1', 'b0' , 'i0', 's0', 'd0' )
        ! We need to handle the variables secondly
        call key_split(key,lkey,rkey)
        if ( len_trim(rkey) == 0 ) then
          lvls = 0
          rkey = lkey
        else
          call lua_open(tbl,lkey,lvls = lvls)
        end if
      end select
      !      print *,'Attempt retrieving: ',trim(key), ' type= ',t
      select case ( t )
      case ( 'a1' )
        call associate(a1,v)
        call lua_get(tbl,rkey,V0)
        na1 = len_trim(V0)
        a1 = ' '
        a1(1:na1) = cpack(V0(1:na1))
      case ( 'b0' ) 
        call associate(b0,v)
        call lua_get(tbl,rkey,b0)
      case ( 'b1' ) 
        call associate(b1,v)
        call lua_get(tbl,key,b1)
      case ( 'b2' ) 
        call associate(b2,v)
        call lua_get(tbl,key,b2)
      case ( 'i0' ) 
        call associate(i0,v)
        call lua_get(tbl,rkey,i0)
      case ( 'i1' ) 
        call associate(i1,v)
        call lua_get(tbl,key,i1)
      case ( 'i2' ) 
        call associate(i2,v)
        call lua_get(tbl,key,i2)
      case ( 's0' ) 
        call associate(s0,v)
        call lua_get(tbl,rkey,s0)
      case ( 's1' ) 
        call associate(s1,v)
        call lua_get(tbl,key,s1)
      case ( 's2' ) 
        call associate(s2,v)
        call lua_get(tbl,key,s2)
      case ( 'd0' ) 
        call associate(d0,v)
        call lua_get(tbl,rkey,d0)
        !         print *,'getting: '//trim(key)//' as ',d0
      case ( 'd1' ) 
        call associate(d1,v)
        call lua_get(tbl,key,d1)
      case ( 'd2' ) 
        call associate(d2,v)
        call lua_get(tbl,key,d2)
        !         print *,'getting: '//trim(key)//' as ',d2
      end select
      !      print *,'Done retrieving: ',trim(key), ' type= ',t
      call lua_close(tbl,lvls = lvls)
    end subroutine get_var

  end subroutine flook_if_get_dict

  subroutine key_split(key,lkey,rkey)
    character(len=*), intent(in) :: key
    character(len=*), intent(inout) :: lkey, rkey
    integer :: i
    i = index(key,'.',back=.true.)
    if ( i > 0 ) then
      lkey = trim(adjustl(key(1:i-1)))
      rkey = key(i+1:)
    else
      lkey = trim(adjustl(key))
      rkey = ' '
    end if
  end subroutine key_split


  function flook_if_print_objects(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use esl_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(dict) :: et
    character(len=DICT_KEY_LENGTH) :: key
    character(len=12) :: fmt = '(tr2,a,'','')'

    ! Currently we only let the current io-node
    ! print out information.
    nret = 0
    if ( .not. flook_if_debug ) return

    ! Print out information
    write(*,'(a)') '-- esl table structure available in LUA'
    write(*,'(a)') 'esl = {'
    write(*,'(tr2,a,i0,'','')') 'Node = ',Node + 1

    ! Loop across all keys in the dictionaries
    et = .first. options
    do while ( .not. (.empty. et) )
      key = .key. et
      write(*,fmt) trim(key)
      et = .next. et
    end do

    ! Loop across all keys in the dictionaries
    et = .first. variables
    do while ( .not. (.empty. et) )
      key = .key. et
      write(*,fmt) trim(key)
      et = .next. et
    end do

    write(*,'(a)') '}'

  end function flook_if_print_objects

#endif

end module esl_flook_if_m
