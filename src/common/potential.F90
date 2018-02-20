module esl_potential_m

  use prec, only : dp,ip
  use esl_energy_m
  use esl_geometry_m
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
    real(dp), allocatable :: vxc(:)  !xc potential

    real(dp) :: ionicOffset !< Offset of the external potential

    type(psolver_t) :: psolver
    type(xc_t) :: xc
  contains
    procedure, public :: init
    procedure, public :: calculate
    procedure, public :: compute_ext_loc
    final  :: cleanup
  end type potential_t

contains

  !Initialize the potentials
  !----------------------------------------------------
  subroutine init(pot, grid, states, geo, periodic)
    class(potential_t) :: pot
    type(grid_t),    intent(in) :: grid
    type(states_t),  intent(in) :: states
    type(geometry_t),intent(in) :: geo
    logical,         intent(in) :: periodic
    
    character(len = 1) :: geocode

    pot%np = grid%np

    allocate(pot%hartree(1:pot%np))
    pot%hartree(1:pot%np) = 0.d0
    allocate(pot%external(1:pot%np))
    pot%external(1:pot%np) = 0.d0
    allocate(pot%vxc(1:pot%np))
    pot%vxc(1:pot%np) = 0.d0

    if (periodic) then
      geocode = 'P'
    else
      geocode = 'F'
    end if

    pot%ionicOffset = 0._dp
    call pot%psolver%init(0, 1, geocode, grid%ndims, grid%hgrid)

    call pot%xc%init(geo, grid)

    call pot%compute_ext_loc(grid, geo)

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
    real(dp),           intent(in)    :: density(:)
    type(energy_t),     intent(inout) :: energy

    integer :: i

    forall (i = 1:pot%np)
      pot%hartree(i) = density(i)
    end forall

    call pot%psolver%hartree_potential(pot%hartree, pot%external, pot%ionicOffset, energy%hartree)

    !Here we need to compute the xc potential
    call pot%xc%calculate(density, energy, pot%vxc)

  end subroutine calculate


  !Compute the local part of the external potential
  !----------------------------------------------------
  subroutine compute_ext_loc(pot, grid, geo)
    use esl_constants_m, only : PI
    class(potential_t), intent(inout) :: pot
    type(grid_t),       intent(in)    :: grid
    type(geometry_t),   intent(in)    :: geo

    integer :: iat, ip, is
    real(kind=dp), allocatable :: extloc(:)

    allocate(extloc(1:pot%np)) 
    pot%external(1:pot%np) = 0.d0
 
    do iat = 1, geo%n_atoms
      is = geo%species_idx(iat)

      ! Convert the radial density to the cartesian grid
      call grid%radial_function(geo%species(is)%vlocal, geo%xyz(:,iat), func=extloc)
      
      do ip=1, pot%np
        pot%external(ip) = pot%external(ip) + extloc(ip)
      end do
    end do

    deallocate(extloc)

  end subroutine compute_ext_loc

end module esl_potential_m
