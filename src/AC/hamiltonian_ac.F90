!< This routine handles the calculation of various parts of the sparse Hamiltonian matrix
!<
!< The input sparse matrix *must* be pre-allocated.
module esl_hamiltonian_ac_m

  use prec, only: dp
  use pspiof_m
  use fdf, only: fdf_get
  
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
  subroutine calculate_H0(this, basis, geom, grid)
    class(hamiltonian_ac_t), intent(inout) :: this
    type(basis_ac_t), intent(in) :: basis
    type(geometry_t), intent(in) :: geom
    type(grid_t), intent(in) :: grid

    ! Calculate individual elements for the H0-elements
    this%kin%M = 0._dp
    call hamiltonian_ac_laplacian(basis, grid, this%kin)

    ! Calculate V_kb matrix elements
    this%vkb%M(:) = 0._dp
    call hamiltonian_ac_Vkb(geom, basis, grid, this%vkb)

  end subroutine calculate_H0

  !< Constructs the initial H from the H0 terms
  !<
  !< This is using pre-calculated values:
  !<   H = H_kin + H_vkb
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
  
  subroutine hamiltonian_ac_Vkb(geom, basis, grid, H)
    class(geometry_t), intent(in) :: geom
    class(basis_ac_t), intent(in) :: basis
    type(grid_t), intent(in) :: grid
    type(sparse_matrix_t), intent(inout) :: H

    ! Local variables
    type(sparse_pattern_t), pointer :: sp
    ! The projector
    type(pspiof_meshfunc_t) :: proj
    ! Grid parts
    real(dp), allocatable :: iG(:), pG(:), jG(:)
    real(dp) :: cut_off

    ! Loop basis
    integer :: ib, ibs, io, iio, il, im
    real(dp) :: ibxyz(3), ir_max

    ! Loop KB
    integer :: a, as, ap, apl, m
    real(dp) :: axyz(3), apr_max
    real(dp), allocatable :: KB_i(:)
    real(dp) :: Ep

    ! Loop j
    integer :: ind, jo, jb, jbs, jjo, jl, jm
    real(dp) :: jbxyz(3), jr_max

    ! Final calculation
    real(dp) :: Vkb

    ! Immediately return, if not needed
    if ( .not. H%initialized() ) return

    ! The option for figuring out KB projector overlap is
    ! determined by the cutoff radius for the projectors.
    ! In this case we limit the projectors to the
    ! cutoff value where the smallest R value where:
    !    r_max @ abs(F(R)) < CutOff
    cut_off = fdf_get('Basis.AC.KB.Cutoff', 0.00001_dp)

    ! Retrieve pointer
    sp => H%sp

    ! Allocate the grid
    allocate(iG(grid%np))
    allocate(pG(grid%np))
    allocate(jG(grid%np))

    ! Loop over all orbital connections in the sparse pattern and
    ! once we have an i,j we loop over atoms to find projectors within a
    ! close range.

    ! loop: ib (basis sites)
    basis_loop: do ib = 1, basis%n_site
      ibs = basis%site_state_idx(ib)
      ibxyz = basis%xyz(:, ib)

      ! loop: i (row in H)
      i_loop: do io = basis%site_orbital_start(ib), basis%site_orbital_start(ib + 1) - 1
        ! Orbital index on basis site
        iio = io - basis%site_orbital_start(ib) + 1

        ! Retrieve the current basis-functions
        !   r_max maximum radius
        !   l quantum number
        !   m quantum number
        ir_max = basis%state(ibs)%orb(iio)%r_cut
        il = basis%state(ibs)%orb(iio)%l
        im = basis%state(ibs)%orb(iio)%m
        
        ! Calculate the basis function on the grid -> iG
        call grid%radial_function_ylm(basis%state(ibs)%orb(iio)%R, il, im, ibxyz, iG)

        ! loop: alpha
        ! Since there are typically few projectors we will loop those first
        ! that should ease the sorting of which projectors are close to ibxyz
        
        ! If 
        !   VKB_ij == 0 since <KB_alpha|phi_i> == 0
        atom_loop: do a = 1, geom%n_atoms
          as = geom%species_idx(a)
          axyz = geom%xyz(:, a)

          ! Loop projector
          KB_loop: do ap = 1, geom%species(as)%n_projectors

            ! Retrieve the r_max
            apr_max = geom%species(as)%get_projector_rmax(ap, cut_off)

            ! Now we check whether ap is non-zero
            if ( not_within_cutoff(ibxyz, ir_max, axyz, apr_max) ) cycle

            ! Retrieve the projector
            call geom%species(as)%get_projector(ap, proj, Ep, apl)

            ! We have a match!
            ! Calculate
            !    <KB_alpha|phi_i> * Ep

            ! Allocate the <KB_alpha|phi_i> for all m quantum numbers of
            ! this KB projector
            allocate(KB_i(-apl:apl))
            do m = -apl, apl

              call grid%radial_function_ylm(proj, apl, m, axyz, pG)

              KB_i(m) = grid%overlap(ibxyz, iG, ir_max, axyz, pG, apr_max) * Ep

            end do

            ! loop: j
            do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1
              
              ! Figure out which basis orbital this belongs too
              jo = sp%column(ind)
              ! Figure out the basis site of the orbital
              jb = basis%orbital_site(jo)
              jbxyz = basis%xyz(:, jb)
              jbs = basis%site_state_idx(jb)
              jjo = jo - basis%site_orbital_start(jb) + 1
              
              ! Retrieve the current basis-functions
              !   r_max maximum radius
              !   l quantum number
              !   m quantum number
              jr_max = basis%state(jbs)%orb(jjo)%r_cut
              jl = basis%state(jbs)%orb(jjo)%l
              jm = basis%state(jbs)%orb(jjo)%m

              ! Now we check whether ap is non-zero
              if ( not_within_cutoff(axyz, apr_max, jbxyz, jr_max) ) cycle

              ! Calculate the basis function on the grid -> jG
              call grid%radial_function_ylm(basis%state(jbs)%orb(jjo)%R, jl, jm, jbxyz, jG)

              Vkb = 0._dp
              do m = -apl, apl
                
                call grid%radial_function_ylm(proj, apl, m, axyz, pG)
                
                Vkb = Vkb + KB_i(m) * grid%overlap(jbxyz, jG, jr_max, axyz, pG, apr_max)
                
              end do
              
              ! Add element to the matrix
              H%M(ind) = H%M(ind) + Vkb

              ! DEBUG print
              ! NOTE this will print out multiple times per diagonal
              !      element corresponding to the number of KB with overlap.
!              if ( ibs == jbs .and. iio == jjo ) &
!                  print *,'# Diagonal Vkb matrix: ', ibs, iio, H%M(ind)
              
            end do
            
            ! Clean-up
            deallocate(KB_i)

            call pspiof_meshfunc_free(proj)
            
          end do KB_LOOP
          
        end do atom_loop
      end do i_loop
      
    end do basis_loop
    
    ! Clean-memory
    deallocate(iG, pG, jG)
    
  contains
    
    pure function not_within_cutoff(xyz1, r1, xyz2, r2) result(not_within)
      real(dp), intent(in) :: xyz1(3), r1
      real(dp), intent(in) :: xyz2(3), r2
      real(dp) :: d
      logical :: not_within

      d = (xyz1(1) - xyz2(1)) ** 2 + &
          (xyz1(2) - xyz2(2)) ** 2 + &
          (xyz1(3) - xyz2(3)) ** 2
      
      not_within = d > (r1 + r2) ** 2
      
    end function not_within_cutoff
    
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
