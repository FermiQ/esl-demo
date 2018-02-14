module esl_hamiltonian_m
  use prec, only : dp,ip
  use esl_basis_m
  use esl_geometry_m
  use esl_grid_m
  use esl_hamiltonian_pw_m
  use esl_potential_m
  use esl_states_m

#ifdef WITH_MPI
  use mpi
#endif


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
    procedure, public :: eigensolver
    final :: cleanup
  end type hamiltonian_t

contains

  !Initialize the Hamiltonian
  !----------------------------------------------------
  subroutine init(this, basis, geo, states, periodic)
    class(hamiltonian_t) :: this
    type(basis_t),    intent(in) :: basis
    type(geometry_t), intent(in) :: geo
    type(states_t),   intent(in) :: states
    logical,          intent(in) :: periodic
 
    integer mpicomm

    call this%energy%init()
    call this%potentials%init(basis%grid, states, geo, periodic)

    mpicomm = 0
#ifdef WITH_MPI
    mpicomm = MPI_COMM_WORLD
#endif

    select case (basis%type)
    case ( PLANEWAVES )
      call this%hm_pw%init(this%potentials, mpicomm)
    case ( ATOMCENTERED )
    !TODO
    end select
 

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
      call this%hm_pw%eigensolver(states, basis%pw)
    case ( ATOMCENTERED )
    !TODO
    end select

  end subroutine eigensolver


  
end module esl_hamiltonian_m
