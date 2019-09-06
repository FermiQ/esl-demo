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
  use esl_states_m

#ifdef WITH_MPI
  use mpi_dist_block_cyclic_m
#endif
  use esl_elsi_m, only: elsi_t
  use esl_calc_density_matrix_ac_m

  implicit none

  private

  public :: density_ac_t

  !Data structure for the density
  type density_ac_t
    !< Density on the real space grid
    real(dp), allocatable :: rho(:)

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
    procedure, public :: calculate_density_matrix
    procedure, nopass, public :: add_density_matrix
    final  :: finalizer

  end type density_ac_t

contains

  !< Initialize the density matrices for this object
  subroutine init(this, basis)
    class(density_ac_t), intent(inout) :: this
    type(basis_ac_t), intent(in) :: basis

    allocate(this%rho(1:basis%grid%np))

    ! Initialize the sparse matrices
    call this%DM%init(basis%sparse_pattern)
    call this%EDM%init(basis%sparse_pattern)

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine finalizer(this)
    type(density_ac_t), intent(inout) :: this

    if( allocated(this%rho) ) then
      deallocate(this%rho)
    end if

    call this%DM%delete()
    call this%EDM%delete()

  end subroutine finalizer

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

    ! We initialize the guess on the grid here
    this%rho(:) = 0._dp
    call add_density_matrix(basis, this%DM, this%rho)

  end subroutine guess

  !< Calculate the density from the hosted density matrix in this object
  subroutine calculate(this, basis)
    !< Density matrix and real-space density
    class(density_ac_t), intent(inout) :: this
    !< Atomic orbital basis
    type(basis_ac_t), intent(in) :: basis

    this%rho(:) = 0._dp
    call add_density_matrix(basis, this%DM, this%rho)

  end subroutine calculate

  !< Calculate a new density matrix from an input states object
  subroutine calculate_density_matrix(this, states)
    !< Density matrix and real-space density
    class(density_ac_t), intent(inout) :: this
    !< Atomic orbital basis
    type(states_t), intent(in) :: states

    type(sparse_pattern_t), pointer :: sp
    integer :: ispin, ikpt, io, ind, jo, ib
    integer :: nr, nc
    real(dp) :: kw
    real(dp) :: occ_k, rio_occ, erio_occ
    complex(dp) :: cio_occ, ecio_occ

    sp => this%DM%sp
    nr = sp%nr
    nc = sp%nc

    this%DM%M(:) = 0._dp
    this%EDM%M(:) = 0._dp

    do ikpt = 1, states%nkpt
      kw = states%k_weights(ikpt)
      do ispin = 1, states%nspin
        
        do ib = 1, states%nstates
          ! Calculate new DM
          if ( states%complex_states ) then
            print *,'density_ac::calculate_DM::complex to be implemented!'
          end if
          occ_k = states%occ_numbers(ib,ispin,ikpt) * kw
          
          do io = 1, nr
            if ( states%complex_states ) then
              cio_occ = states%states(ib,ispin,ikpt)%zcoef(io) * occ_k
              ecio_occ = cio_occ * states%eigenvalues(ib,ispin,ikpt)
            else
              rio_occ = states%states(ib,ispin,ikpt)%dcoef(io) * occ_k
              erio_occ = rio_occ * states%eigenvalues(ib,ispin,ikpt)
            end if

            ! Loop columns
            do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1
              jo = sp%column(ind)
              if ( states%complex_states ) then
                ! TODO phases from k
                !this%DM%M(ind) = this%DM%M(ind) + real(cio_occ * dconjg(states%states(ib,ispin,ikpt)%zcoef(jo)), dp)
                !this%EDM%M(ind) = this%EDM%M(ind) + real(ecio_occ * dconjg(states%states(ib,ispin,ikpt)%zcoef(jo)), dp)
              else
                this%DM%M(ind) = this%DM%M(ind) + rio_occ * states%states(ib,ispin,ikpt)%dcoef(jo)
                this%EDM%M(ind) = this%EDM%M(ind) + erio_occ * states%states(ib,ispin,ikpt)%dcoef(jo)
              end if

            end do
            
          end do
        end do
      end do
    end do
    
  end subroutine calculate_density_matrix

  !< Calculate and print-out the Mulliken charges
  subroutine mulliken_summary(this, basis, S)
    use yaml_output

    type(density_ac_t), intent(in) :: this
    type(basis_ac_t), intent(in) :: basis
    type(sparse_matrix_t), intent(in) :: S

    ! Local variables for calculating the Mulliken charges
    integer :: ia, io, ind
    integer :: io1, io2

    type(sparse_pattern_t), pointer :: sp
    
    ! Actual Mulliken charges per site
    real(dp), allocatable :: M(:), F(:)
    character(len=16) :: str

    ! Retrieve sparse pattern
    sp => S%sp
    
    allocate(M(basis%n_site))
    allocate(F(basis%n_orbital))

    ! Calculate the total number charge, per site
    do ia = 1 , basis%n_site
      
      M(ia) = 0._dp
      do io = basis%site_orbital_start(ia), basis%site_orbital_start(ia + 1) - 1

        F(io) = 0._dp
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

          F(io) = F(io) + S%M(ind) * this%DM%M(ind)
          
        end do

        M(ia) = M(ia) + F(io)

      end do

    end do

    ! Produce YAML output
    call yaml_mapping_open("Mulliken")
    call yaml_map('Total', sum(M))
    call yaml_comment("Q", hfill = "-")
    do ia = 1, basis%n_site

      write(str, '(i0)') ia
      call yaml_mapping_open(trim(str))
      !str = 'StillNotCreated' !basis%species(basis%site_state_idx(ia))%label
      !call yaml_map('Label', trim(str))
      call yaml_map('Sum', M(ia))
      io1 = basis%site_orbital_start(ia)
      io2 = basis%site_orbital_start(ia + 1) - 1
      call yaml_map('Individual', F(io1:io2))
      call yaml_mapping_close()
      
    end do

    call yaml_mapping_close()

    deallocate(M, F)
    
  end subroutine mulliken_summary
    

  !< Add a sparse density matrix to the density grid using basis coefficients, etc.
  !<
  !< Add the basis functions density to the grid via an input density matrix.
  !TODO change to allow for spin-polarized calculations (DM(2), rho(:,2))
  subroutine add_density_matrix(basis, DM, rho)

    !< Grid container that defines this density object
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

    loop_grid: do ip = 1, basis%grid%np

      ! Initialize the density to be added to this grid-point
      rho_ip = 0._dp

      ! Figure out all neighbouring sites
      call populate_neighbours(n_neigh, basis%grid%r(:, ip))

      ! Loop on all sites connecting to the grid-point
      loop_neigh: do in = 1, n_neigh

        ! Retrieve current site
        ia = neigh(in)
        is = basis%site_state_idx(ia)

        ! Retrieve current vector from grid-point to the basis-orbital
        idr(:) = basis%grid%r(:,ip) - basis%xyz(:,ia)

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

            jdr(:) = basis%grid%r(:,ip) - basis%xyz(:,ja)
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
