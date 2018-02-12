!< Module for handling density matrices
!<
!< This module should implement the basics of creating
!< an initial (valence filled) density matrix and/or
!< read-in/extrapolate a DM to a different structure.
module esl_density_matrix_ac_m

  use prec
  use esl_basis_ac_m, only: basis_ac_t
  use esl_grid_m, only: grid_t
  use esl_sparse_pattern_m, only: sparse_pattern_t
  use esl_sparse_matrix_m, only: sparse_matrix_t

  implicit none

  public :: density_matrix_ac_next
  public :: density_matrix_ac_init_atomic
  public :: density_matrix_ac_grid_charge

contains

  !< Prepare the next density matrix
  !<
  !< This routine gets passed an old sparse pattern and a
  !< new sparse pattern.
  !< When the old sparse pattern is not allocated (created)
  !< we automatically initialize the DM with the atomic fillings.
  subroutine density_matrix_ac_next(basis, old_sp, new_sp, DM)

    class(basis_ac_t), intent(in) :: basis
    type(sparse_pattern_t), intent(in) :: old_sp
    type(sparse_pattern_t), intent(in), target :: new_sp
    ! One element per spin component.
    type(sparse_matrix_t), intent(inout), allocatable :: DM(:)

    if ( old_sp%initialized() ) then

      ! For now we still do the atomic fillings...
      call density_matrix_ac_init_atomic(basis, new_sp, DM)

    else

      call density_matrix_ac_init_atomic(basis, new_sp, DM)

    end if

  end subroutine density_matrix_ac_next

  !< Initialize the diagonal density matrix with atomic fillings
  subroutine density_matrix_ac_init_atomic(basis, sp, DM)

    use prec, only: dp

    class(basis_ac_t), intent(in) :: basis
    type(sparse_pattern_t), intent(in), target :: sp
    ! One element per spin component.
    type(sparse_matrix_t), intent(inout), allocatable :: DM(:)

    real(dp) :: frac_s
    integer :: nspin
    integer :: ia, is, io, iio, ind

    nspin = 1
    ! Initialize the density matrix
    if ( allocated(DM) ) then
      ! We assume it is allocated correctly (nspin is correct)
    else
      ! TODO insert correct size  (number of spin-components)
      allocate(DM(nspin))
    end if

    do is = 1, nspin
      call DM(is)%init(sp)

      ! Initialize all elements to zero
      DM(is)%M = 0._dp
    end do

    ! Define the fraction of filling per spin channel
    frac_s = 1._dp / nspin

    ! Loop over all orbital connections in the sparse pattern and
    ! set the diagonal density matrix
    do ia = 1, basis%n_sites
      is = basis%species_idx(ia)

      ! Loop on orbitals
      do io = basis%site_function_start(ia), basis%site_function_start(ia + 1) - 1
        ! Orbital index on atom
        iio = io - basis%site_function_start(ia) + 1

        ! Loop entries in the sparse pattern
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1
          ! Skip if it does not correspond to the diagonal element
          if ( sp%column(ind) /= io ) cycle

          ! Set diagonal component
          ! TODO add atomic charge specification
          !          DM(:)%M(ind) = basis% TODO basis %info(is)%q0(iio) * frac_s

        end do

      end do

    end do

  end subroutine density_matrix_ac_init_atomic


  !< Add the basis functions density to the grid via the density matrix
  subroutine density_matrix_ac_grid_charge(basis, DM, grid, rho)

    use esl_numeric_m, only: grylmr

    !< Atomic orbital basis
    class(basis_ac_t), intent(in) :: basis
    !< Density matrix
    type(sparse_matrix_t), intent(in) :: DM
    !< Grid to expand the density matrix on
    class(grid_t), intent(inout) :: grid
    !< Associated density on the grid
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

      ! Figure out all neighbouring sites
      call populate_neighbours(n_neigh, grid%r(:, ip))

      ! Loop on all sites connecting to the grid-point
      loop_neigh: do in = 1, n_neigh

        ! Retrieve current site
        ia = neigh(in)
        is = basis%species_idx(ia)

        ! Retrieve current vector from grid-point to the basis-orbital
        idr(:) = grid%r(:,ip) - basis%sites_xyz(:,ia)
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

            jdr(:) = grid%r(:,ip) - basis%sites_xyz(:,ja)
            ! TODO retrieve quantum numbers l and m
            call grylmr(jdr(1), jdr(2), jdr(3), jl, jm, jao)

            ! Add density contribution to the grid
            rho(ip) = rho(ip) + iao * DM%M(ind) * jao

          end do

        end do loop_DM

      end do loop_neigh

    end do loop_grid

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

        dist = sqrt(sum( (basis%sites_xyz(:, is) - r) ** 2 ))

        ! TODO add check for r-max
        if ( dist <= r_max ) then
          n_neigh = n_neigh + 1
          neigh(n_neigh) = is
        end if

      end do

    end subroutine populate_neighbours

  end subroutine density_matrix_ac_grid_charge

end module esl_density_matrix_ac_m
