!< Dictionary container to enable ALL dictionaries containing variables
module esl_dict_m

  use prec, only : dp
#ifdef WITH_FLOOK
  use dictionary
#endif

  implicit none

#ifdef WITH_FLOOK

  private

  ! A dictionary for all variables
  type(dict), public :: esl_variables

  interface esl_dict_var_add
    module procedure dict_variable_add_v_0d
!    module procedure dict_variable_add_a_1d
    module procedure dict_variable_add_b_0d
    module procedure dict_variable_add_i_0d
    module procedure dict_variable_add_i_1d
    module procedure dict_variable_add_d_0d
    module procedure dict_variable_add_d_1d, dict_variable_add_d_2d
  end interface esl_dict_var_add

  ! Public interface
  public :: esl_dict_var_add

contains

  subroutine dict_variable_add_v_0d(name,val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.true.)
    esl_variables = esl_variables // (name.kv.trim(val))
  end subroutine dict_variable_add_v_0d
!  subroutine dict_variable_add_a_1d(name,val)
!    character(len=*), intent(in) :: name
!    character(len=1), intent(inout), target :: val(:)
!    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
!    esl_variables = esl_variables // (name.kvp.val)
!  end subroutine dict_variable_add_a_1d
  subroutine dict_variable_add_b_0d(name,val)
    character(len=*), intent(in) :: name
    logical, intent(inout), target :: val
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
    esl_variables = esl_variables // (name.kvp.val)
  end subroutine dict_variable_add_b_0d
  subroutine dict_variable_add_i_0d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout), target :: val
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
    esl_variables = esl_variables // (name.kvp.val)
  end subroutine dict_variable_add_i_0d
  subroutine dict_variable_add_i_1d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout), target :: val(:)
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
    esl_variables = esl_variables // (name.kvp.val)
  end subroutine dict_variable_add_i_1d
  subroutine dict_variable_add_d_0d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
    esl_variables = esl_variables // (name.kvp.val)
  end subroutine dict_variable_add_d_0d
  subroutine dict_variable_add_d_1d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val(:)
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
    esl_variables = esl_variables // (name.kvp.val)
  end subroutine dict_variable_add_d_1d
  subroutine dict_variable_add_d_2d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout), target :: val(:,:)
    if ( name.in.esl_variables ) call delete(esl_variables,name,dealloc=.false.)
    esl_variables = esl_variables // (name.kvp.val)
  end subroutine dict_variable_add_d_2d

#endif

end module esl_dict_m

