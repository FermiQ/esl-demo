module esl_density_ac_m

  use prec, only : dp
  use esl_basis_ac_m
  use esl_energy_m
  use esl_grid_m
  use esl_mulliken_ac_m
  use esl_potential_m
  use esl_sparse_pattern_m
  use esl_sparse_matrix_m
  use esl_hamiltonian_ac_m

  use mpi_dist_block_cyclic_m
  use esl_elsi_m, only: elsi_t
  use esl_calc_density_matrix_ac_m
  

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
    this%EDM%M(:) = 0._dp

  end subroutine guess

  !< Calculate the density from the hosted density matrix in this object
  subroutine calculate(this, elsi, grid, pot, basis, S, rho, energy, out)
#ifdef WITH_MPI
    use mpi, only: mpi_comm_world
!    include 'mpif.h'
#endif
    class(density_ac_t), intent(inout) :: this
    !< ELSI handler
    class(elsi_t), intent(inout) :: elsi
    !< Grid container that defines this density object
    class(grid_t), intent(in) :: grid
    !< Potential container which calculates Hartree and XC
    class(potential_t), intent(inout) :: pot
    !< Atomic orbital basis
    class(basis_ac_t), intent(in) :: basis
    !< The overlap matrix (has to be pre-calculated on entry)
    type(sparse_matrix_t) :: S
    !< Real-space grid on which to calculate the (total) density from the DM
    real(dp), intent(inout) :: rho(:)
    !< Energy type to contain all energies
    type(energy_t), intent(inout) :: energy
    !< Output density
    type(density_ac_t), intent(inout) :: out

    !< Real-space grid on which to retain the atomic filled density
    real(dp), allocatable :: rho_atom(:)

    !< Create the initial DM
    type(sparse_matrix_t) :: DM_atom

    ! Currently we contain the Hamiltonian here.
    ! However, since it is needed elsewhere we should decide where to place it.
    type(sparse_matrix_t) :: H
    type(mpi_dist_block_cyclic_t) :: dist

    ! Initialize the distribution
    call dist%init(MPI_COMM_World, this%DM%sp%nr, this%DM%sp%nr)

    ! TODO logic for calculating the output density from an input
    ! density.
    allocate(rho_atom(grid%np))
    ! Initialize
    rho_atom(:) = 0._dp
    
    ! Calculate the atomic density
!    call basis%atomic_density_matrix(DM_atom)
    
    ! Expand the atomic dM on an auxiliary grid
!    call add_density_matrix(grid, basis, DM_atom, rho_atom)
!    call DM_atom%delete()

!    print *, 'DEBUG rho-atom sum', grid%integrate(rho_atom)

    ! 1. Start by calculating the density from the DM on the grid
    rho(:) = 0._dp
    call add_density_matrix(grid, basis, this%DM, rho)

!    rho_atom = rho - rho_atom
!    print *, 'DEBUG dRho sum', grid%integrate(rho_atom)

    ! Initialize the Hamiltonian to 0
    call H%init(S%sp)
    H%M(:) = 0._dp

    ! First step for any SCF step is to setup the
    ! correct Hamiltonian

    !  1. Add the Laplacian to the Hamiltonian
    call hamiltonian_ac_laplacian(basis, grid, H)
    ! Calculate the kinetic energy
    energy%kinetic = sum(H%M * this%DM%M)

    !  2. Add Hartree potential
    ! Now we are in a position to calculate Hartree potential from
    ! rho

    print *, 'DEBUG Rho sum', grid%integrate(rho)

    ! Now call the potentials_t%calculate which does:
    !   1. Calculate the XC potential.
    !   2. Calculate the Hartree potential
    !   3. Calculate the external local potential
    call pot%calculate(rho, energy)

    ! Sum external, Hartree and XC potential,
    ! this will ease the addition of the matrix elements to do it in one go.
    ! Re-use rho-atom as the sum of potentials.
    ! I.e. after this line we cannot use rho_atom anymore!
    rho_atom(:) = pot%hartree(:) + pot%vxc(:) + pot%external(:)
    print *, 'DEBUG V sum', grid%integrate(rho_atom)
    call hamiltonian_ac_potential(basis, grid, rho_atom, H)

    ! Calculate the eigenvalue energy
    energy%eigenvalues = sum(H%M * this%DM%M)

    ! X. Calculate output density matrix elements from the Hamiltonian
!    call set_elsi_sparsity_pattern_ac(elsi, dist, H%sp)
!    call calc_density_matrix_ac(elsi, H, S, out%DM, energy%fermi)

    call dist%delete()

    call my_check()

    ! Final, output Mulliken charges
    call mulliken_ac_summary(basis, S, out%DM)
    call energy%display()

  contains

    subroutine my_check()

      real(dp), allocatable :: eig(:)
      real(dp), allocatable :: B(:), work(:)
      integer :: n, info
      integer :: io, jo, ib, ind

      type(sparse_pattern_t), pointer :: sp

      integer :: N_Ef

      sp => H%sp
      n = sp%nr

      allocate(eig(n))
      allocate(B(n ** 2))
      allocate(work(n ** 2))
      B = S%M
      
      call dsygv(1, 'V', 'U', n, H%M, n, B, n, eig, work, n ** 2,info)

      N_Ef = nint( basis%Q )
      energy%fermi = (eig(n_ef) + eig(n_ef+1)) / 2
      print *, 'DEBUG fermi level: ', energy%fermi, eig

      ! Re-construct the DM
      out%DM%M(:) = 0._dp

      ! Loop rows
      do io = 1, n
        
        ! Loop columns
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1
          jo = sp%column(ind)
          
          ! Loop bands
          do ib = 1, N_ef
            
            ! Calculate new DM
            out%DM%M(ind) = out%DM%M(ind) + &
                H%M((ib-1)*n + jo) * H%M((ib-1)*n + io)
            
          end do
          
        end do
        
      end do
      
      deallocate(B,eig,work)

    end subroutine my_check

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
