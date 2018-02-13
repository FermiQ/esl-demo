module esl_density_pw_m

  use prec, only : dp
  use esl_basis_pw_m
  use esl_geometry_m
  use esl_grid_m
  use esl_message_m
  use esl_states_m
  use esl_utils_pw_m
  use yaml_output

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
    procedure, public :: residue
    procedure, public :: get_den
    procedure, public :: set_den
    final  :: cleanup
  end type density_pw_t

contains

  !Initialize the density
  !----------------------------------------------------
  subroutine init(this, grid, basis)
    class(density_pw_t),     intent(inout) :: this
    type(grid_t),               intent(in) :: grid
    type(basis_pw_t), target,   intent(in) :: basis

    if(grid%np /= product(basis%ndims)) then
      call message_error("Number of grid point and number of plane waves are not consistent.")
    end if

    allocate(this%density(1:grid%np))
    this%np = grid%np

    this%pw => basis

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, geo, grid)
    class(density_pw_t), intent(inout) :: this
    type(geometry_t), intent(in) :: geo
    type(grid_t),     intent(in) :: grid
    
    real(dp), allocatable :: atomicden(:)
    integer :: iat, ip, is
    real(dp):: norm

    this%density(1:grid%np) = 0.d0

    allocate(atomicden(1:grid%np))
    atomicden(1:grid%np) = 0.d0

    call yaml_mapping_open("Guess atomic density")

    ! We expect only atoms to contain initial density
    do iat = 1, geo%n_atoms
      is = geo%species_idx(iat)      

      !Convert the radial density to the cartesian grid
      call grid%radial_function(geo%species(is)%rho, 0, 0, geo%xyz(:,iat), atomicden)

      call integrate(grid, atomicden, norm)
      call yaml_map("Norm", norm)

      !Summing up to the total density
      forall (ip = 1:grid%np)
        this%density(ip) = this%density(ip) + atomicden(ip)
      end forall
    end do

    deallocate(atomicden)

    call yaml_mapping_close()

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
      kpt(1:3) = 0.d0
      do isp = 1, states%nspin
        do ist = 1, states%nstates
          !TODO: Here we should have a gmap for each k-point
          !From the G vectors to the real space
          call pw2grid(this%pw%grid, this%pw%gmap, this%pw%ndims, this%pw%npw, &
                               states%states(ik,isp,ist)%zcoef, coef_rs)
 
          !We accumulate the density
          do ip = 1, this%np
            this%density(ip) = this%density(ip) + real(coef_rs(ip)*conjg(coef_rs(ip)),kind=dp)
          end do
        end do
      end do
    end do

    deallocate(coef_rs)

  end subroutine calculate

  !Residue
  !---------------------------------------------------
  function residue(this, grid, nel, other) result(res)
    class(density_pw_t), intent(in) :: this
    type(grid_t),        intent(in) :: grid
    integer,             intent(in) :: nel
    type(density_pw_t),  intent(in) :: other

    real(kind=dp), allocatable :: diff(:)
    real(kind=dp) :: res
    integer :: ip

    allocate(diff(1:grid%np)) 

    !We use rhonew to compute the relative density
    do ip = 1, this%np
      diff(ip) = abs(other%density(ip) - this%density(ip))
    end do
    call integrate(grid, diff, res)
    res = res/real(nel)

    deallocate(diff)

  end function residue

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
