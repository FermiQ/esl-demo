module esl_density_t

  use prec, only : dp

  use esl_basis_m
  use esl_grid_m
  use esl_system_m

  implicit none
  private

  public :: density_t

  !Data structure for the density
  type density_t
     integer :: np !< Copied from grid

     real(dp), allocatable :: density(:)
     real(dp), allocatable :: density_matrix(:,:)
   contains

     procedure, public :: init
     procedure, public :: guess
     procedure, public :: calculate
     procedure, public :: get_den
     procedure, public :: set_den
     final  :: cleanup

  end type density_t

contains

  !Initialize the density
  !----------------------------------------------------
  subroutine init(this, basis, grid)
    class(density_t), intent(inout) :: this
    type(basis_t), intent(in) :: basis
    type(grid_t), intent(in) :: grid

    !Parse the informations from the input file
    select case(basis%type)
    case(PLANEWAVES)
       !Initialization structures for the PW case
       allocate(this%density(1:grid%np))
       this%density(1:grid%np) = 0.d0
    case(ATOMICORBS)
       !Initialization structures for the LO case
       !TEMP
       allocate(this%density(1:grid%np))
       this%density(1:grid%np) = 0.d0
    end select

    this%np = grid%np

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, system)
    class(density_t), intent(inout) :: this
    type(system_t), intent(in) :: system

    real(dp), allocatable :: radial(:)
    real(dp), allocatable :: atomicden(:)
    integer :: iat, np_radial, ip

    this%density(1:system%grid%np) = 0.d0

    allocate(atomicden(1:system%grid%np))
    atomicden(1:system%grid%np) = 0.d0

    do iat = 1, system%natoms
       !Get the number of points in the radial grid
       np_radial = 1
       allocate(radial(1:np_radial))
       !Get atomic density on the radial grid

       !Convert the radial density to the cartesian grid

       !We do not need the radial density anymore
       deallocate(radial)

       !Summing up to the total density
       forall(ip=1:system%grid%np)
          this%density(ip) = this%density(ip) + atomicden(ip)
       end forall
    end do

    deallocate(atomicden)

  end subroutine guess

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(density_t), intent(inout) :: this

    if(allocated(this%density_matrix)) deallocate(this%density_matrix)
    if(allocated(this%density)) deallocate(this%density)

  end subroutine cleanup

  !Calc density
  !----------------------------------------------------
  subroutine calculate(this, basis)
    class(density_t), intent(inout) :: this
    type(basis_t), intent(in) :: basis

    select case(basis%type)
    case(PLANEWAVES)

    case(ATOMICORBS)

    end select

  end subroutine calculate


  !Copy the density to an array
  !----------------------------------------------------
  subroutine get_den(this, rho)
    class(density_t) :: this
    real(dp), intent(out) :: rho(:)

    integer :: ip

    forall(ip=1:this%np)
       rho(ip) = this%density(ip)
    end forall

  end subroutine get_den

  !Copyi the density from an array
  !----------------------------------------------------
  subroutine set_den(this, rho)
    class(density_t) :: this
    real(dp), intent(in) :: rho(:)

    integer :: ip

    forall(ip=1:this%np)
       this%density(ip) = rho(ip)
    end forall

  end subroutine set_den

end module esl_density_t
