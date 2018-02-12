module esl_density_pw_m

  use prec, only : dp
  use esl_basis_pw_m
  use esl_geometry_m
  use esl_grid_m
  use esl_message_m
  use esl_states_m
  use esl_utils_pw_m

  implicit none

  private

  public :: density_pw_t

  !Data structure for the density
  type density_pw_t
    integer :: np !< Copied from grid

    real(dp), allocatable :: density(:)
    type(basis_pw_t), pointer :: pw
  contains
    private
    procedure, public :: init
    procedure, public :: guess
    procedure, public :: calculate
    procedure, public :: get_den
    procedure, public :: set_den
    final  :: cleanup
  end type density_pw_t

contains

  !Initialize the density
  !----------------------------------------------------
  subroutine init(this, grid, basis_pw)
    class(density_pw_t),     intent(inout) :: this
    type(grid_t),               intent(in) :: grid
    type(basis_pw_t), target,   intent(in) :: basis_pw

    if(grid%np /= product(basis_pw%ndims)) then
      call message_error("Number of grid point and number of plane waves are not consistent.")
    end if

    allocate(this%density(1:grid%np))
    this%density(1:grid%np) = 0.d0

    this%np = grid%np

    this%pw => basis_pw

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, geo, grid)
    class(density_pw_t), intent(inout) :: this
    type(geometry_t), intent(in) :: geo
    type(grid_t),     intent(in) :: grid
    
    real(dp), allocatable :: radial(:)
    real(dp), allocatable :: atomicden(:)
    integer :: iat, np_radial, ip

    this%density(1:grid%np) = 0.d0

    allocate(atomicden(1:grid%np))
    atomicden(1:grid%np) = 0.d0

    ! We expect only atoms to contain initial density
    do iat = 1, geo%n_atoms
      
      !Get the number of points in the radial grid
      np_radial = 1
      allocate(radial(1:np_radial))
      !Get atomic density on the radial grid

      !Convert the radial density to the cartesian grid

      !We do not need the radial density anymore
      deallocate(radial)

      !Summing up to the total density
      forall (ip = 1:grid%np)
        this%density(ip) = this%density(ip) + atomicden(ip)
      end forall
    end do

    deallocate(atomicden)

  end subroutine guess

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(density_pw_t), intent(inout) :: this

    if(allocated(this%density)) deallocate(this%density)
    nullify(this%pw)

  end subroutine cleanup

  !Calc density
  !----------------------------------------------------
  subroutine calculate(this, states)
    class(density_pw_t), intent(inout) :: this
    type(states_t),   intent(in) :: states
    
    integer :: ik, isp, ist, ip
    complex(kind=dp), allocatable :: coef_rs(:)
 
    real(kind=dp) :: kpt(3)

    !Coefficient in real space
    allocate(coef_rs(1:this%np))

    this%density(:) = 0.d0

    ! Density should be calculated from states
    do ik = 1, states%nkpt
     kpt(3) = 0.d0
      do isp = 1, states%nspin
        do ist = 1, states%nstates
          !From the G vectors to the real space
          call fourier2grid(this%pw%grid, this%pw%gmet, kpt, this%pw%ndims, this%pw%ecut, &
                               states%states(ik,isp,ist)%zcoef, this%np, coef_rs)
 
          !We accumulate the density
          do ip = 1, this%np
            this%density(ip) = this%density(ip) + real(coef_rs(ip)*conjg(coef_rs(ip)),kind=dp)
          end do
        end do
      end do
    end do

    deallocate(coef_rs)

  end subroutine calculate


  !Copy the density to an array
  !----------------------------------------------------
  subroutine get_den(this, rho)
    class(density_pw_t) :: this
    real(dp), intent(out) :: rho(:)

    integer :: ip

    forall(ip = 1:this%np)
      rho(ip) = this%density(ip)
    end forall

  end subroutine get_den

  !Copyi the density from an array
  !----------------------------------------------------
  subroutine set_den(this, rho)
    class(density_pw_t) :: this
    real(dp), intent(in) :: rho(:)

    integer :: ip

    forall (ip = 1:this%np)
      this%density(ip) = rho(ip)
    end forall

  end subroutine set_den

end module esl_density_pw_m
