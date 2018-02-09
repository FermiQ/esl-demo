module potential_esl
  use prec, only : dp,ip

  use basis_esl
  use density_esl
  use energy_esl
  use grid_esl
  use psolver_esl
  use states_esl
  use xc_esl

  implicit none
  private

  public :: potential_t

  !Data structure for the potentials
  type potential_t
     integer :: np !< Number of points in real space

     real(kind=dp), allocatable :: hartree(:)  !Hartree potential
     real(kind=dp), allocatable :: external(:) !External local potential
     real(kind=dp), allocatable :: vxc(:,:)  !xc potential

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
  subroutine init(pot, basis, grid, states)
    class(potential_t) :: pot
    type(basis_t), intent(in) :: basis
    type(grid_t),  intent(in) :: grid
    type(states_t), intent(in):: states

    character(len = 1) :: geocode

    pot%np = grid%np

    allocate(pot%hartree(1:pot%np))
    allocate(pot%external(1:pot%np))
    allocate(pot%vxc(1:pot%np, 1:states%nspin))

    select case(basis%basis_type)
    case(PLANEWAVES)
       geocode = 'P'
    case(ATOMICORBS)
       geocode = 'F'
    end select

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
    type(density_t),   intent(in)    :: density
    type(energy_t),    intent(inout) :: energy

    call pot%psolver%h_potential(density, pot%hartree, pot%np, &
         & pot%external, pot%ionicOffset, energy%hartree)

    !Here we need to compute the xc potential
    call pot%xc%calculate(density, pot%vxc)

  end subroutine calculate

end module potential_esl
