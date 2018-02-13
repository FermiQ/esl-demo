module esl_density_ac_m

  use prec, only : dp
  use esl_basis_ac_m
  use esl_grid_m
  use esl_sparse_pattern_m
  use esl_sparse_matrix_m

  implicit none

  private

  public :: density_ac_t

  !Data structure for the density
  type density_ac_t

    !< Density matrix
    type(sparse_matrix_t) :: DM
    
    !< Energy density matrix
    type(sparse_matrix_t) :: EDM

  contains

    private
    procedure, public :: init
    procedure, public :: residue

    procedure, public :: guess
    procedure, public :: calculate
    procedure, public :: add_density_matrix
    final  :: cleanup

  end type density_ac_t

contains

  !Initialize the density. Allocate the full grid
  !----------------------------------------------------
  subroutine init(this)
    class(density_ac_t), intent(inout) :: this

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(density_ac_t), intent(inout) :: this

    call this%DM%delete()
    call this%EDM%delete()

  end subroutine cleanup

  !< Guess the initial density from the atomic orbitals
  !<
  !< The guessed fillings are based on the basis%atomic_density_matrix
  !< which is the neutral atom valence fillings.
  subroutine guess(this, basis)
    class(density_ac_t), intent(inout) :: this
    !< Atomic orbital basis
    type(basis_ac_t), intent(in) :: basis

    ! Create the atomic density matrix
    call basis%atomic_density_matrix(this%DM)

  end subroutine guess

  !< Calculate the density from the hosted density matrix in this object
  subroutine calculate(this, grid, basis, rho, out)
    class(density_ac_t), intent(inout) :: this
    !< Grid container that defines this density object
    class(grid_t), intent(in) :: grid
    !< Atomic orbital basis
    class(basis_ac_t), intent(in) :: basis
    !< Real-space grid on which to calculate the density from the DM
    real(dp), intent(inout) :: rho(:)
    !< Output density
    type(density_ac_t), intent(inout) :: out

    ! TODO logic for calculating the output density from an input
    ! density.

    ! I.e. we should calculate the Hamiltonian etc.
    
    ! This should contain
    call this%add_density_matrix(grid, basis, this%DM, rho)
    
  end subroutine calculate
    

  !< Add a sparse density matrix to the density grid using basis coefficients, etc.
  !<
  !< Add the basis functions density to the grid via an input density matrix.
  subroutine add_density_matrix(this, grid, basis, DM, rho)

    use esl_numeric_m, only: grylmr

    class(density_ac_t), intent(inout) :: this
    !< Grid container that defines this density object
    class(grid_t), intent(in) :: grid
    !< Atomic orbital basis
    class(basis_ac_t), intent(in) :: basis
    !< Density matrix
    type(sparse_matrix_t), intent(in) :: DM
    !< Real-space density
    real(dp), intent(inout) :: rho(:)

    ! Local variables
    type(sparse_pattern_t), pointer :: sp
    integer :: ip
    !< Neighbour list that contains a consecutive list of species
    !< neighbouring a grid point
    integer, allocatable :: neigh(:)
    integer :: n_neigh ! current number of neighbours

    ! Loop counters
    integer :: in, ind
    ! First site indices
    integer :: ia, is, io, iio
    ! Second site indices
    integer :: ja, js, jo, jjo

    ! Looping species quantum numbers (l and m)
    integer :: il, im, jl, jm
    ! Looping species distance to the grid-point and the atomic orbital at the grid-point
    real(dp) :: idr(3), iao, jdr(3), jao

    ! Total density at grid-point
    real(dp) :: rho_ip

    ! Retrieve the sparse pattern
    sp => DM%sp

    ! The following algorithm will loop over all grid points and
    ! subsequently on all atoms to figure out which atoms are within the range of
    ! their orbital radius.
    ! The routine does not impose periodic boundary conditions but that could be
    ! added by using MOD operations etc.

    ! Initialize the neighbour list
    allocate(neigh(basis%n_sites))

    loop_grid: do ip = 1, grid%np

      ! Initialize the density to be added to this grid-point
      rho_ip = 0._dp

      ! Figure out all neighbouring sites
      call populate_neighbours(n_neigh, grid%r(:, ip))

      ! Loop on all sites connecting to the grid-point
      loop_neigh: do in = 1, n_neigh

        ! Retrieve current site
        ia = neigh(in)
        is = basis%species_idx(ia)

        ! Retrieve current vector from grid-point to the basis-orbital
        idr(:) = grid%r(:,ip) - basis%xyz(:,ia)

        ! TODO retrieve quantum numbers l and m
        call grylmr(idr(1), idr(2), idr(3), il, im, iao)

        ! Loop over basis DM functions
        loop_DM: do io = basis%site_function_start(ia), basis%site_function_start(ia + 1) - 1

          iio = io - basis%site_function_start(ia) + 1 ! local basis-function index of site ia

          ! Loop over density matrix for the neighbouring basis-sites
          ! TODO this loop could be reduced by limiting the basis functions
          ! with respect to individual ranges.
          do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

            ! Figure out which atom this orbital belongs too
            jo = sp%column(ind) ! global basis-function index
            ! Figure out the atomic index of the orbital
            ja = basis%function_site(jo) ! global basis site
            js = basis%species_idx(ja) ! which basis specie
            jjo = jo - basis%site_function_start(ja) + 1 ! local basis-function index of site ja

            jdr(:) = grid%r(:,ip) - basis%xyz(:,ja)
            ! TODO retrieve quantum numbers l and m
            call grylmr(jdr(1), jdr(2), jdr(3), jl, jm, jao)

            ! Add density contribution to the grid
            rho_ip = rho_ip + iao * DM%M(ind) * jao

          end do

        end do loop_DM

      end do loop_neigh

      rho(ip) = rho(ip) + rho_ip

    end do loop_grid

    ! Clean up
    deallocate(neigh)

  contains

    subroutine populate_neighbours(n_neigh, r)
      integer, intent(inout) :: n_neigh
      real(dp), intent(in) :: r(3)
      real(dp) :: dist, r_max
      integer :: is

      n_neigh = 0
      r_max = 10._dp ! TODO retrieve correct radius of basis-functions
      do is = 1, basis%n_sites

        dist = sqrt(sum( (basis%xyz(:, is) - r) ** 2 ))

        ! TODO add check for r-max
        if ( dist <= r_max ) then
          n_neigh = n_neigh + 1
          neigh(n_neigh) = is
        end if

      end do

    end subroutine populate_neighbours

  end subroutine add_density_matrix

  !< Calculate the residual of two density matrices
  function residue(this, other) result(res)
    class(density_ac_t), intent(in) :: this, other

    real(dp) :: res
    
    res = maxval( abs( this%DM%M - other%DM%M ) )

  end function residue

end module esl_density_ac_m
