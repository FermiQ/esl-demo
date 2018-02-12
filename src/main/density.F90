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

    real(dp), allocatable :: rhoin(:)
    real(dp), allocatable :: rhoout(:)
    real(dp), allocatable :: rhonew(:)

    type(mixing_t)   :: mixer
    type(density_ac_t) :: density_ac
    type(density_pw_t) :: density_pw
  contains
    private
    procedure, public :: init
    procedure, public :: guess
    procedure, public :: calculate
    procedure, public :: get_relden
    procedure, public :: mix
    final  :: cleanup
  end type density_t

contains

  !Initialize the density
  !----------------------------------------------------
  subroutine init(this, grid, basis)
    class(density_t), intent(inout) :: this
    type(grid_t),     intent(in) :: grid
    type(basis_t),    intent(in) :: basis

    allocate(this%rhoin(1:grid%np))
    this%rhoin(1:grid%np) = 0.d0

    allocate(this%rhoout(1:grid%np))
    this%rhoin(1:grid%np) = 0.d0

    allocate(this%rhonew(1:grid%np))
    this%rhoin(1:grid%np) = 0.d0

    this%np = grid%np

    call this%mixer%init()

    select case (basis%type)
    case ( PLANEWAVES )
      call this%density_pw%init(grid, basis%pw)
    case ( ATOMCENTERED )
      call this%density_ac%init(grid)
    end select

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, basis, geo, grid)
    class(density_t), intent(inout) :: this
    type(basis_t),    intent(in) :: basis
    type(geometry_t), intent(in) :: geo
    type(grid_t),     intent(in) :: grid
   
    select case (basis%type)
    case ( PLANEWAVES )
      call this%density_pw%guess(geo, grid)
    case ( ATOMCENTERED )
      call this%density_ac%guess(geo, grid)
    end select
 

  end subroutine guess

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(density_t), intent(inout) :: this

    if(allocated(this%rhoin)) deallocate(this%rhoin)
    if(allocated(this%rhoout)) deallocate(this%rhoout)
    if(allocated(this%rhonew)) deallocate(this%rhonew)

  end subroutine cleanup

  !Calc density
  !----------------------------------------------------
  subroutine calculate(this, basis, states)
    class(density_t), intent(inout) :: this
    type(basis_t),       intent(in) :: basis
    type(states_t),      intent(in) :: states

    select case (basis%type)
    case ( PLANEWAVES )
      !Saving the in density for the mixing
      call this%density_pw%get_den(this%rhoin)
      !Calc. density
      call this%density_pw%calculate(states)
      !Saving the out density for the mixing
      call this%density_pw%get_den(this%rhoout)
    case ( ATOMCENTERED )
      !Saving the in density for the mixing
      call this%density_ac%get_den(this%rhoin)
      !Calc. density
      call this%density_ac%calculate()
      !Saving the out density for the mixing
      call this%density_ac%get_den(this%rhoout)
    end select

  end subroutine calculate


  !Copy the relative density
  !----------------------------------------------------
  real(kind=dp) function get_relden(this, grid, nel) result(reldens)
    class(density_t) :: this
    type(grid_t), intent(in) :: grid
    integer,      intent(in) :: nel

    integer :: ip

    !Test tolerance and print status
    !We use rhonew to compute the relative density
    do ip = 1, this%np
      this%rhonew(ip) = abs(this%rhoout(ip) - this%rhoin(ip))
    end do
    call integrate(grid, this%rhonew, reldens)
    reldens = reldens/real(nel)


  end function get_relden

  !Copyi the density from an array
  !----------------------------------------------------
  subroutine mix(this, basis)
    class(density_t) :: this
    type(basis_t),  intent(in) :: basis

    select case (basis%type)
    case ( PLANEWAVES )  
      call mixing_linear(this%mixer, this%np, this%rhoin, this%rhoout, this%rhonew)
      call this%density_pw%set_den(this%rhonew)
    case ( ATOMCENTERED )
      call mixing_linear(this%mixer, this%np, this%rhoin, this%rhoout, this%rhonew)
      call this%density_ac%set_den(this%rhonew)
    end select

  end subroutine mix

end module esl_density_m
