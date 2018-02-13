module esl_density_m

  use prec, only : dp
  use esl_basis_m
  use esl_density_ac_m
  use esl_density_pw_m
  use esl_geometry_m
  use esl_grid_m
  use esl_mixing_m
  use esl_states_m

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

  !Initialize the density
  !----------------------------------------------------
  subroutine init(this, grid, basis)
    class(density_t), intent(inout) :: this
    type(grid_t),     intent(in) :: grid
    type(basis_t),    intent(in) :: basis

    allocate(this%rho(1:grid%np))
    this%rho(1:grid%np) = 0.d0
    this%np = grid%np

    select case ( basis%type )
    case ( PLANEWAVES )
      call this%density_pw%init(grid, basis%pw)
    case ( ATOMCENTERED )
      
      call this%ac%init()
      
    end select

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, basis, geo)
    class(density_t), intent(inout) :: this
    type(basis_t),    intent(in) :: basis
    type(geometry_t), intent(in) :: geo
   
    select case ( basis%type )
    case ( PLANEWAVES )

      ! TODO fix interface 
      call this%density_pw%guess(geo, basis%grid)
      
    case ( ATOMCENTERED )

      ! Guess the initial density from the atomic fillings
      call this%ac%guess(basis%grid, basis%ac)
      
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
  subroutine calculate(this, basis, states, out)
    class(density_t),   intent(inout) :: this
    type(basis_t),       intent(in) :: basis
    type(states_t),      intent(in) :: states
    type(density_t),  intent(inout) ::  out

    ! Calculate density on the grid
    select case (basis%type)
    case ( PLANEWAVES )

<<<<<<< HEAD
      call this%density_pw%calculate(states)
    case ( ATOMCENTERED )

      ! Calculate the density on the grid
      call this%ac%calculate(basis%grid, basis%ac, this%rho, out%ac)
      
    end select

  end subroutine calculate


  !< Calculate the relative difference between two densities
  function residue(this, basis, other) result(res)
    class(density_t), intent(in) :: this, other
    type(basis_t), intent(in) :: basis

    real(dp) :: res

        ! Calculate density on the grid
    select case ( basis%type )
    case ( PLANEWAVES )

      ! TODO implement the residue function in density_pw_t
!      res = this%density_pw%residue(basis%grid, basis%pw, other%density_pw)
      
    case ( ATOMCENTERED )

      res = this%ac%residue(other%ac)

    end select

=======
      ! Calculate density
      call this%density_pw%calculate(states)

    case ( ATOMCENTERED )

      ! Calculate the density on the grid
      call this%ac%calculate(basis%grid, basis%ac, this%rho, out%ac)
      
    end select

  end subroutine calculate


  !< Calculate the relative difference between two densities
  function residue(this, basis, other) result(res)
    class(density_t), intent(in) :: this, other
    type(basis_t), intent(in) :: basis

    real(dp) :: res

        ! Calculate density on the grid
    select case ( basis%type )
    case ( PLANEWAVES )

      ! TODO implement the residue function in density_pw_t
!      res = this%density_pw%residue(basis%grid, basis%pw, other%density_pw)
      
    case ( ATOMCENTERED )

      res = this%ac%residue(other%ac)

    end select

>>>>>>> 4c2f0e6aa9d4f98348842c29c1c0f203ed0a6a55
    
!!$    !Test tolerance and print status
!!$    !We use rhonew to compute the relative density
!!$    do ip = 1, this%np
!!$      this%rhonew(ip) = abs(this%rhoout(ip) - this%rhoin(ip))
!!$    end do
!!$    call integrate(grid, this%rhonew, reldens)
!!$    reldens = reldens/real(nel)

  end function residue

end module esl_density_m
