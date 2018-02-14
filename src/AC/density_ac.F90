module esl_density_ac_m

  use prec, only : dp
  use esl_basis_ac_m
  use esl_grid_m
  use esl_mulliken_ac_m
  use esl_sparse_pattern_m
  use esl_sparse_matrix_m
  use esl_hamiltonian_ac_m

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
    procedure, nopass, public :: add_density_matrix
    final  :: cleanup

  end type density_ac_t

contains

  !< Initialize the density matrices for this object
  subroutine init(this, sp)
    class(density_ac_t), intent(inout) :: this
    type(sparse_pattern_t), target, intent(in) :: sp

    ! Initialize the sparse matrices
    call this%DM%init(sp)
    call this%EDM%init(sp)

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
  subroutine calculate(this, grid, basis, S, rho, out)
    class(density_ac_t), intent(inout) :: this
    !< Grid container that defines this density object
    class(grid_t), intent(in) :: grid
    !< Atomic orbital basis
    class(basis_ac_t), intent(in) :: basis
    !< The overlap matrix (has to be pre-calculated on entry)
    type(sparse_matrix_t) :: S
    !< Real-space grid on which to calculate the (total) density from the DM
    real(dp), intent(inout) :: rho(:)
    !< Output density
    type(density_ac_t), intent(inout) :: out

    !< Real-space grid on which to retain the atomic filled density
    real(dp), allocatable :: rho_atom(:)

    !< Create the initial DM
    type(sparse_matrix_t) :: DM_atom

    ! Currently we contain the Hamiltonian here.
    ! However, since it is needed elsewhere we should decide where to place it.
    type(sparse_matrix_t) :: H

    ! TODO logic for calculating the output density from an input
    ! density.
    allocate(rho_atom(grid%np))
    ! Initialize
    rho_atom(:) = 0._dp
    
    ! Calculate the atomic density
    call basis%atomic_density_matrix(DM_atom)
    
    ! Expand the atomic dM on an auxiliary grid
    call add_density_matrix(grid, basis, DM_atom, rho_atom)
    call DM_atom%delete()

    print *, 'DEBUG rho-atom sum', grid%integrate(rho_atom)

    ! 1. Start by calculating the density from the DM on the grid
    call add_density_matrix(grid, basis, this%DM, rho)

    ! Now begin by calculating the Laplacian for the Hamiltonian
    call H%init(S%sp)
    
    ! Initialize the Hamiltonian to 0
    H%M(:) = 0._dp

    ! Add the Laplacian to the Hamiltonian
    call hamiltonian_ac_laplacian(basis, grid, H)

    ! Now we are in a position to calculate Hartree potential from
    ! rho and/or rho_atom
    

    ! 2. Calculate matrix elements for the Hamiltonian

    call mulliken_ac_summary(basis, S, this%DM)

  end subroutine calculate
    

  !< Add a sparse density matrix to the density grid using basis coefficients, etc.
  !<
  !< Add the basis functions density to the grid via an input density matrix.
  subroutine add_density_matrix(grid, basis, DM, rho)

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
    integer :: max_no ! number of orbitals per state

    ! Looping species quantum numbers (l and m)
    integer :: il, im, jl, jm
    ! Looping species distance to the grid-point and the atomic orbital at the grid-point
    real(dp) :: idr(3), jdr(3), jpsi
    real(dp), allocatable :: ipsi(:)

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
    allocate(neigh(basis%n_site))

    max_no = 0
    do is = 1 , basis%n_state
      max_no = max(max_no, basis%state(is)%n_orbital)
    end do
    
    ! Allocate one array to calculate all psi for the LHS 
    allocate(ipsi(max_no))

    loop_grid: do ip = 1, grid%np

      ! Initialize the density to be added to this grid-point
      rho_ip = 0._dp

      ! Figure out all neighbouring sites
      call populate_neighbours(n_neigh, grid%r(:, ip))

      ! Loop on all sites connecting to the grid-point
      loop_neigh: do in = 1, n_neigh

        ! Retrieve current site
        ia = neigh(in)
        is = basis%site_state_idx(ia)

        ! Retrieve current vector from grid-point to the basis-orbital
        idr(:) = grid%r(:,ip) - basis%xyz(:,ia)

        ! Calculate all psi-components at point r for all basis-functions
        ! on this specie
        call basis%get_psi(is, idr, ipsi)

        ! Loop over basis DM functions
        loop_DM: do io = basis%site_orbital_start(ia), basis%site_orbital_start(ia + 1) - 1

          iio = io - basis%site_orbital_start(ia) + 1 ! local basis-function index of site ia

          ! Loop over density matrix for the neighbouring basis-sites
          
          ! TODO this loop could be reduced by limiting the basis functions
          ! with respect to individual ranges.
          do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

            ! Figure out which atom this orbital belongs too
            jo = sp%column(ind) ! global basis-function index
            ! Figure out the atomic index of the orbital
            ja = basis%orbital_site(jo) ! global basis site
            js = basis%site_state_idx(ja) ! which basis specie
            jjo = jo - basis%site_orbital_start(ja) + 1 ! local basis-function index of site ja

            jdr(:) = grid%r(:,ip) - basis%xyz(:,ja)
            call basis%get_psi(js, jjo, jdr, jpsi)

            ! Add density contribution to the grid
            rho_ip = rho_ip + ipsi(iio) * DM%M(ind) * jpsi

          end do

        end do loop_DM

      end do loop_neigh

      rho(ip) = rho(ip) + rho_ip

    end do loop_grid

    ! Clean up
    deallocate(neigh, ipsi)

  contains

    subroutine populate_neighbours(n_neigh, r)
      integer, intent(inout) :: n_neigh
      real(dp), intent(in) :: r(3)
      real(dp) :: dist, r_cut
      integer :: is

      n_neigh = 0
      do is = 1, basis%n_site

        dist = sqrt(sum( (basis%xyz(:, is) - r) ** 2 ))

        ! Retrieve cut-off
        r_cut = basis%state(basis%site_state_idx(is))%r_cut
        
        if ( dist <= r_cut ) then
          
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
