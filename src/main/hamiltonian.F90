module esl_hamiltonian_m
  use prec, only : dp,ip
  use esl_basis_m
  use esl_energy_m
  use esl_force_m
  use esl_geometry_m
  use esl_grid_m
  use esl_hamiltonian_pw_m
  use esl_ion_interaction_m
  use esl_potential_m
  use esl_states_m

  implicit none
  private

  public ::                          &
      hamiltonian_t

  !Data structure for the Hamiltonian
  type hamiltonian_t
    type(energy_t)          :: energy
    type(force_t)           :: force
    type(ion_interaction_t) :: ion_inter
    type(potential_t)       :: potentials
    type(hamiltonian_pw_t)  :: hm_pw
  contains
    private
    procedure, public :: init
    procedure, public :: eigensolver
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

    call this%energy%init()
    call this%potentials%init(grid, states, geo, periodic)
    call this%force%init(geo%n_atoms)
    call this%ion_inter%init()

    if (periodic) then
      call this%ion_inter%calculate_periodic(geo, this%force%ionion, this%energy%ionion)
    else
      call this%ion_inter%calculate_isolated(geo, this%force%ionion, this%energy%ionion)
    end if

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(hamiltonian_t) :: this

  end subroutine cleanup

  !Eigensolver
  !----------------------------------------------------
  subroutine eigensolver(this, basis, states)
    class(hamiltonian_t) :: this
    type(basis_t),  intent(in)    :: basis
    type(states_t), intent(inout) :: states

    select case (basis%type)
    case ( PLANEWAVES )
      call this%hm_pw%eigensolver(states)
    case ( ATOMCENTERED )
    !TODO
    end select

  end subroutine eigensolver


  
end module esl_hamiltonian_m
