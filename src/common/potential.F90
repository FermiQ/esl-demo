module esl_potential_m

  use prec, only : dp,ip
  use esl_energy_m
  use esl_grid_m
  use esl_psolver_m
  use esl_states_m
  use esl_xc_m

  implicit none
  private

  public :: potential_t

  !Data structure for the potentials
  type potential_t
    integer :: np !< Number of points in real space

    real(dp), allocatable :: hartree(:)  !Hartree potential
    real(dp), allocatable :: external(:) !External local potential
    real(dp), allocatable :: vxc(:,:)  !xc potential

    real(dp) :: ionicOffset !< Offset of the external potential

    type(psolver_t) :: psolver
    type(xc_t) :: xc
  contains
    procedure, public :: init
    procedure, public :: calculate
    final  :: cleanup
  end type potential_t

contains

  !Initialize the potentials
  !----------------------------------------------------
  subroutine init(pot, grid, states, periodic)
    class(potential_t) :: pot
    type(grid_t),   intent(in) :: grid
    type(states_t), intent(in) :: states
    logical,        intent(in) :: periodic
    
    character(len = 1) :: geocode

    pot%np = grid%np

    allocate(pot%hartree(1:pot%np))
    allocate(pot%external(1:pot%np))
    allocate(pot%vxc(1:pot%np, 1:states%nspin))

    if (periodic) then
      geocode = 'P'
    else
      geocode = 'F'
    end if

    pot%ionicOffset = 0._dp
    call pot%psolver%init(0, 1, geocode, grid%ndims, grid%hgrid)

    call pot%xc%init()
    !Here we need to init the libxc and pspio parts

  end subroutine init


  !Release the potentials
  !----------------------------------------------------
  subroutine cleanup(pot)
    type(potential_t):: pot

    if(allocated(pot%hartree)) deallocate(pot%hartree)
    if(allocated(pot%external)) deallocate(pot%external)
    if(allocated(pot%vxc)) deallocate(pot%vxc)

  end subroutine cleanup


  !Compute the different potentials
  !----------------------------------------------------
  subroutine calculate(pot, density, energy)
    class(potential_t), intent(inout) :: pot
    real(dp),           intent(in)    :: density(:,:)
    type(energy_t),     intent(inout) :: energy

    integer :: i

    forall (i = 1:pot%np)
      pot%hartree(i) = sum(density(i,:))
    end forall

    call pot%psolver%h_potential(pot%hartree, pot%external, pot%ionicOffset, energy%hartree)

    !Here we need to compute the xc potential
    call pot%xc%calculate(density, pot%vxc)

  end subroutine calculate

end module esl_potential_m
