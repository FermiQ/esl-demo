module esl_hamiltonian_m
  use prec, only : dp,ip
  use esl_basis_m
  use esl_geometry_m
  use esl_grid_m
  use esl_hamiltonian_pw_m
  use esl_potential_m
  use esl_states_m

  implicit none
  private

  public ::                          &
      hamiltonian_t

  !Data structure for the Hamiltonian
  type hamiltonian_t
    type(potential_t)       :: potential
    type(hamiltonian_pw_t)  :: hm_pw
  contains
    private
    procedure, public :: init
    final :: cleanup
  end type hamiltonian_t

contains

  !Initialize the Hamiltonian
  !----------------------------------------------------
  subroutine init(this, grid, geo, states, periodic)
    class(hamiltonian_t) :: this
    type(grid_t),     intent(in) :: grid
    type(geometry_t), intent(in) :: geo
    type(states_t),   intent(in) :: states
    logical,          intent(in) :: periodic

    call this%potential%init(grid, states, geo, periodic)

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(hamiltonian_t) :: this

  end subroutine cleanup

  
end module esl_hamiltonian_m
