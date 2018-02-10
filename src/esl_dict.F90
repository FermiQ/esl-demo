!< Dictionary container to enable ALL dictionaries containing variables
module esl_dict_m

  use prec, only : dp
#ifdef WITH_LUA
  use dictionary
#endif

  implicit none

#ifdef WITH_LUA

  private

  ! A dictionary for all the options
  type(dict), public :: esl_options

  ! A dictionary for all variables
  type(dict), public :: esl_variables

  interface esl_dict_add
    module procedure dict_variable_add_v_0d
    module procedure dict_variable_add_a_1d
    module procedure dict_variable_add_b_0d
    module procedure dict_variable_add_i_0d
    module procedure dict_variable_add_i_1d
    module procedure dict_variable_add_d_0d
    module procedure dict_variable_add_d_1d, dict_variable_add_d_2d
  end interface esl_dict_add

  ! Public interface
  public :: esl_dict_clean
  public :: esl_dict_populate, esl_dict_populate_options, esl_dict_populate_variables
  public :: esl_dict_add

contains

  subroutine esl_dict_clean()
    
    call delete(esl_options, dealloc=.false.)
    call delete(esl_variables, dealloc=.false.)
    
  end subroutine esl_dict_clean
  
  subroutine esl_dict_populate()
    
    call esl_dict_populate_options()
    call esl_dict_populate_variables()
    
  end subroutine esl_dict_populate

  subroutine esl_dict_populate_options()

    ! We simply re-create the options, (note the 
    ! de-allocation by "nullification")
    call delete(esl_options, dealloc=.false.)
    
  end subroutine esl_dict_populate_options

  subroutine esl_dict_populate_variables()
    
    ! We simply re-create the options, (note the 
    ! de-allocation by "nullification")
    call delete(esl_variables, dealloc=.false.)

  end subroutine esl_dict_populate_variables

  subroutine dict_variable_add_v_0d(name,val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.true.)
    variables = variables // (name.kv.trim(val))
  end subroutine dict_variable_add_v_0d
  subroutine dict_variable_add_a_1d(name,val)
    character(len=*), intent(in) :: name
    character(len=1), intent(inout), target :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_a_1d
  subroutine dict_variable_add_b_0d(name,val)
    character(len=*), intent(in) :: name
    logical, intent(inout), target :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_b_0d
  subroutine dict_variable_add_i_0d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout), target :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_i_0d
  subroutine dict_variable_add_i_1d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout), target :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_i_1d
  subroutine dict_variable_add_d_0d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_0d
  subroutine dict_variable_add_d_1d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_1d
  subroutine dict_variable_add_d_2d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val(:,:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_2d

#endif

end module esl_dict_m

