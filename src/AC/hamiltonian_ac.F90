!< This routine handles the calculation of various parts of the sparse Hamiltonian matrix
!<
!< The input sparse matrix *must* be pre-allocated.
module esl_hamiltonian_ac_m

  use prec, only: dp
  use esl_basis_ac_m
  use esl_geometry_m
  use esl_grid_m

  use esl_sparse_matrix_m, only: sparse_matrix_t
  use esl_sparse_pattern_m, only: sparse_pattern_t

  implicit none
  
  private

  public :: hamiltonian_ac_laplacian
  public :: hamiltonian_ac_potential
  
  public :: hamiltonian_ac_t

  !< Data structure for the sparse matrix Hamiltonians
  type hamiltonian_ac_t

    ! Store Hamiltonian quantities
    
    type(sparse_matrix_t) :: kin !< Kinetic (Laplacian) Hamiltonian
    type(sparse_matrix_t) :: vkb !< Kleynman-Bylander projectors part of the Hamiltonian
    type(sparse_matrix_t) :: SCF !< SCF Hamiltonian

  contains
    
    procedure, public :: init
    procedure, public :: calculate_H0
    procedure, public :: setup_H0
    
    final :: cleanup
  end type hamiltonian_ac_t

contains

  !< Initialize the Hamiltonian by allocating initial quantities
  !<
  !< This routine will pre-allocate 3 matrices:
  !<  1. The kinetic Hamiltonian (non-SCF dependent)
  !<  2. Non-local Hamiltonian (non-SCF dependent)
  !<  3. SCF Hamiltonian (changed on every SCF cycle)
  subroutine init(this, sparse_pattern)
    class(hamiltonian_ac_t), intent(inout) :: this
    type(sparse_pattern_t), intent(in), target :: sparse_pattern
    
    call this%kin%init(sparse_pattern)
    call this%vkb%init(sparse_pattern)
    call this%SCF%init(sparse_pattern)
    
  end subroutine init
  
  subroutine cleanup(this)
    type(hamiltonian_ac_t), intent(inout) :: this
    
    call this%kin%delete()
    call this%vkb%delete()
    call this%SCF%delete()
    
  end subroutine cleanup

  !< Calculate all non-SCF dependent Hamiltonian terms
  !<
  !< These includes:
  !<  1. The kinetic Hamiltonian (non-SCF dependent)
  !<  2. Non-local Hamiltonian (non-SCF dependent)
  subroutine calculate_H0(this, basis, geo, grid)
    class(hamiltonian_ac_t), intent(inout) :: this
    type(basis_ac_t), intent(in) :: basis
    type(geometry_t), intent(in) :: geo
    type(grid_t), intent(in) :: grid

    ! Calculate individual elements for the H0-elements
    this%kin%M = 0._dp
    call hamiltonian_ac_laplacian(basis, grid, this%kin)

    ! Calculate V_kb matrix elements
    this%vkb%M = 0._dp
    call hamiltonian_ac_Vkb(geo, grid, this%vkb)

  end subroutine calculate_H0

  !< Constructs the initial H from the H0 terms
  !<
  !< This is basically performing:
  !<   H = H_kin + H_vk
  subroutine setup_H0(this)
    class(hamiltonian_ac_t), intent(inout) :: this

    this%SCF%M(:) = this%kin%M(:) + this%vkb%M(:)
    
  end subroutine setup_H0

  subroutine hamiltonian_ac_laplacian(basis, grid, H)
    class(basis_ac_t), intent(in) :: basis
    type(grid_t), intent(in) :: grid
    type(sparse_matrix_t), intent(inout) :: H

    integer :: ia, is, io, iio, ind, jo, ja, js, jjo
    real(dp), allocatable :: iT(:,:), jT(:,:)
    real(dp) :: ixyz(3), ir_max, jxyz(3), jr_max
    integer :: il, im, jl, jm

    type(sparse_pattern_t), pointer :: sp

    ! Immediately return, if not needed
    if ( .not. H%initialized() ) return

    sp => H%sp

    ! Allocate the Laplacian matrices
    allocate(iT(3,grid%np))
    allocate(jT(3,grid%np))

    ! Loop over all orbital connections in the sparse pattern and
    ! calculate the overlap matrix for each of them
    do ia = 1, basis%n_site
      is = basis%site_state_idx(ia)
      ixyz = basis%xyz(:, ia)

      ! Loop on orbitals
      do io = basis%site_orbital_start(ia), basis%site_orbital_start(ia + 1) - 1
        ! Orbital index on atom
        iio = io - basis%site_orbital_start(ia) + 1

        ir_max = basis%state(is)%orb(iio)%r_cut
        il = basis%state(is)%orb(iio)%l
        im = basis%state(is)%orb(iio)%m
        call grid%radial_function_ylm_gradient(basis%state(is)%orb(iio)%R, il, im, ixyz(:), iT)

        ! Loop entries in the sparse pattern
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

          ! Figure out which atom this orbital belongs too
          jo = sp%column(ind)
          ! Figure out the atomic index of the orbital
          ja = basis%orbital_site(jo)
          jxyz = basis%xyz(:, ja)
          js = basis%site_state_idx(ja)
          jjo = jo - basis%site_orbital_start(ja) + 1

          ! We are now in a position to calculate the
          ! overlap matrix. I.e. we know the atom, the
          ! orbital indices and their positions
          jr_max = basis%state(js)%orb(jjo)%r_cut
          jl = basis%state(js)%orb(jjo)%l
          jm = basis%state(js)%orb(jjo)%m
          call grid%radial_function_ylm_gradient(basis%state(js)%orb(jjo)%R, jl, jm, jxyz(:), jT)

          H%M(ind) = H%M(ind) + matrix_T(iT, jT)

          ! DEBUG print
!          if ( ia == ja .and. iio == jjo ) &
!              print *,'# Diagonal kinetic matrix: ', ia, iio, H%M(ind)

        end do

      end do

    end do

    deallocate(iT, jT)

  contains

    pure function matrix_T(iT, jT) result(T)
      use esl_constants_m, only: PI
      real(dp), intent(in) :: iT(:,:), jT(:,:)
      real(dp) :: T
      integer :: ip

      T = 0._dp
      do ip = 1, grid%np
        T = T + iT(1,ip) * jT(1,ip) + iT(2,ip) * jT(2,ip) + iT(3,ip) * jT(3,ip)
      end do
      ! TODO check Laplacian and units, it isn't fully correct, but closer to Siesta (in its current state)
      T = T * grid%volelem / (2*PI)
      
    end function matrix_T

  end subroutine hamiltonian_ac_laplacian
  
  subroutine hamiltonian_ac_Vkb(geo, grid, H)
    class(geometry_t), intent(in) :: geo
    type(grid_t), intent(in) :: grid
    type(sparse_matrix_t), intent(inout) :: H

    type(sparse_pattern_t), pointer :: sp

    ! Immediately return, if not needed
    if ( .not. H%initialized() ) return

    ! Retrieve pointer
    sp => H%sp

  end subroutine hamiltonian_ac_Vkb


  subroutine hamiltonian_ac_potential(basis, grid, pot, H)
    use prec, only: dp
    use esl_basis_ac_m, only: basis_ac_t
    use esl_grid_m, only: grid_t
    use esl_sparse_pattern_m, only: sparse_pattern_t
    use esl_sparse_matrix_m, only: sparse_matrix_t

    !< AC basis used
    class(basis_ac_t), intent(in) :: basis
    !< Current grid information
    type(grid_t), intent(in) :: grid
    !< Potential of which to add the matrix elements to the Hamiltonian
    real(dp), intent(in) :: pot(:)
    !< Hamiltonian to add the matrix elements too
    type(sparse_matrix_t), intent(inout) :: H

    integer :: ia, is, io, iio, ind, jo, ja, js, jjo
    real(dp), allocatable :: ipsi(:), jpsi(:)
    real(dp) :: ixyz(3), ir_max, jxyz(3), jr_max
    integer :: il, im, jl, jm

    type(sparse_pattern_t), pointer :: sp

    ! Immediately return, if not needed
    if ( .not. H%initialized() ) return

    sp => H%sp

    ! Allocate the Laplacian matrices
    allocate(ipsi(grid%np))
    allocate(jpsi(grid%np))

    ! Loop over all orbital connections in the sparse pattern and
    ! calculate the overlap matrix for each of them
    do ia = 1, basis%n_site
      is = basis%site_state_idx(ia)
      ixyz = basis%xyz(:, ia)

      ! Loop on orbitals
      do io = basis%site_orbital_start(ia), basis%site_orbital_start(ia + 1) - 1
        ! Orbital index on atom
        iio = io - basis%site_orbital_start(ia) + 1

        ir_max = basis%state(is)%orb(iio)%r_cut
        il = basis%state(is)%orb(iio)%l
        im = basis%state(is)%orb(iio)%m
        call grid%radial_function_ylm(basis%state(is)%orb(iio)%R, il, im, ixyz, ipsi)

        ! Loop entries in the sparse pattern
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

          ! Figure out which atom this orbital belongs too
          jo = sp%column(ind)
          ! Figure out the atomic index of the orbital
          ja = basis%orbital_site(jo)
          jxyz = basis%xyz(:, ja)
          js = basis%site_state_idx(ja)
          jjo = jo - basis%site_orbital_start(ja) + 1

          ! We are now in a position to calculate the
          ! overlap matrix. I.e. we know the atom, the
          ! orbital indices and their positions
          jr_max = basis%state(js)%orb(jjo)%r_cut
          jl = basis%state(js)%orb(jjo)%l
          jm = basis%state(js)%orb(jjo)%m
          call grid%radial_function_ylm(basis%state(js)%orb(jjo)%R, jl, jm, jxyz, jpsi)

          H%M(ind) = H%M(ind) + &
              grid%matrix_elem(ixyz, ipsi, ir_max, pot, jxyz, jpsi, jr_max)
          
          ! DEBUG print
!          if ( ia == ja .and. iio == jjo ) &
!              print *,'# Diagonal potential matrix: ', ia, iio, H%M(ind)

        end do

      end do

    end do

    deallocate(ipsi, jpsi)

  end subroutine hamiltonian_ac_potential

end module esl_hamiltonian_ac_m
