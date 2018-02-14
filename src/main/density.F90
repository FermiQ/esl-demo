module esl_density_m

  use prec, only : dp
  use esl_basis_m
  use esl_density_ac_m
  use esl_density_pw_m
  use esl_geometry_m
  use esl_mixing_m
  use esl_states_m
  use esl_system_m

  implicit none

  private

  public :: density_t

  !Data structure for the density
  type density_t
    integer :: np !< Copied from grid

    real(dp), allocatable :: rho(:)

    type(density_ac_t) :: ac
    type(density_pw_t) :: density_pw
    
  contains
    
    private
    procedure, public :: init
    procedure, public :: guess
    procedure, public :: calculate
    procedure, public :: residue
    final  :: cleanup
    
  end type density_t

contains

  !< Initialize density objects
  subroutine init(this, system)
    class(density_t), intent(inout) :: this
    type(system_t), intent(in) :: system

    allocate(this%rho(1:system%basis%grid%np))
    this%rho(1:system%basis%grid%np) = 0.d0
    this%np = system%basis%grid%np

    select case ( system%basis%type )
    case ( PLANEWAVES )
      
      call this%density_pw%init(system%basis%grid, system%basis%pw)
      
    case ( ATOMCENTERED )
      
      call this%ac%init(system%sparse_pattern)
      
    end select

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, system)
    class(density_t), intent(inout) :: this
    type(system_t), intent(in) :: system
   
    select case ( system%basis%type )
    case ( PLANEWAVES )

      ! TODO fix interface 
      call this%density_pw%guess(system%geo, system%basis%grid)
      
    case ( ATOMCENTERED )

      ! Guess the initial density from the atomic fillings
      call this%ac%guess(system%basis%ac)
      
    end select

  end subroutine guess

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(density_t), intent(inout) :: this

    if( allocated(this%rho) ) then
      deallocate(this%rho)
    end if

  end subroutine cleanup

  ! Calculate output density from an input density
  !----------------------------------------------------
  subroutine calculate(this, system, states, out)
    class(density_t),   intent(inout) :: this
    type(system_t), intent(inout) :: system
    type(states_t),      intent(in) :: states
    type(density_t),  intent(inout) ::  out

    ! Calculate density on the grid
    select case ( system%basis%type )
    case ( PLANEWAVES )

      ! Calculate density
      call out%density_pw%calculate(states)

    case ( ATOMCENTERED )

      ! Calculate the density on the grid
      call this%ac%calculate(system%basis%grid, system%basis%ac, system%S, this%rho, &
          system%energy, out%ac)
      
    end select

  end subroutine calculate


  !< Calculate the relative difference between two densities
  function residue(this, basis, other, states) result(res)
    class(density_t),    intent(in) :: this, other
    type(basis_t),       intent(in) :: basis
    type(states_t),      intent(in) :: states

    real(dp) :: res

        ! Calculate density on the grid
    select case ( basis%type )
    case ( PLANEWAVES )

      res = this%density_pw%residue(basis%grid, states%nel, other%density_pw)
      
    case ( ATOMCENTERED )

      res = this%ac%residue(other%ac)

    end select

  end function residue

end module esl_density_m
