!< Module to initialize a sparse pattern by inspection of the atomic orbitals
!<
!< The sparse pattern will be determined by the radii of two atomic orbitals.
!< If the distance between the centers of the atomic orbitals are within R_1 + R_2
!< we have a sparse element and will add a sparse element.
!<
!< Note that this will *only* create the initial sparse pattern and will finish
!< by finalizing the sparse_pattern_t
module esl_create_sparse_pattern_ac_m

  implicit none

  public :: create_sparse_pattern_ac_create

contains

  !< The arguments of this routine are
  !< @param system the system with containing atomic coordinates and their basis functions
  !< @param sp the sparse pattern which we need to create, it will be deleted.
  subroutine create_sparse_pattern_ac_create(basis, sp)

    use prec, only: dp
    use esl_basis_ac_m, only: basis_ac_t
    use esl_sparse_pattern_m, only: sparse_pattern_t

    class(basis_ac_t), intent(in) :: basis
    type(sparse_pattern_t), intent(inout) :: sp
    integer :: no, max_no, io, jo
    integer :: ia, ja, is, js
    real(dp) :: ir_cut, jr_cut, r_cut, dist

    ! Start by deallocation of the sparse pattern
    call sp%delete()

    ! Figure out the number of orbitals
    max_no = 0
    do is = 1 , basis%n_state
      max_no = max(max_no, basis%state(is)%n_orbital)
    end do
    ! Total number of basis-functions
    no = basis%n_orbital

    ! Now re-initialize the sparse matrix.
    ! In this case we will assume a maximum of 20 atomic connections
    call sp%init(no, no, np=max_no * 20)

    ! Loop over all sites
    do ia = 1, basis%n_site
      is = basis%site_state_idx(ia)
      ir_cut = basis%state(is)%r_cut

      ! Add the connections to it-self
      call add_elements(ia, ia, 0._dp)

      ! Only loop the remaining atoms (no need to double process)
      do ja = ia + 1, basis%n_site
        js = basis%site_state_idx(ja)
        jr_cut = basis%state(js)%r_cut

        ! Calculate whether the distance between the two
        ! atoms is within their basis range.
        r_cut = ir_cut + jr_cut

        ! Calculate the distance between the two atomic centers
        dist = sqrt( sum((basis%xyz(:,ia) - basis%xyz(:,ja)) ** 2) )
        
        ! Only process if the maximum distance is within range.
        if ( dist <= r_cut ) then

          ! Add all orbitals to the sparse pattern
          call add_elements(ia, ja, dist)
          call add_elements(ja, ia, dist)

        end if

      end do
    end do

    ! Reduce sparse pattern to contained elements
    call sp%finalize()

  contains

    subroutine add_elements(ia, ja, dist)
      integer, intent(in) :: ia, ja
      real(dp), intent(in) :: dist
      real(dp) :: ir_cut, jr_cut
      integer :: isite, jsite
      integer :: is, io, js, jo

      ! TODO, only do it one way (it is symmetric by definition, these floating point operations may!)

      ! Loop orbitals on both atoms
      do io = basis%site_orbital_start(ia) , basis%site_orbital_start(ia+1) - 1
        
        isite = basis%orbital_site(io)
        is = basis%site_state_idx(isite)
        ir_cut = basis%state(is)%orb(io - basis%site_orbital_start(isite) + 1)%r_cut
        
        do jo = basis%site_orbital_start(ja) , basis%site_orbital_start(ja+1) - 1
          
          jsite = basis%orbital_site(jo)
          js = basis%site_state_idx(jsite)
          jr_cut = basis%state(js)%orb(jo - basis%site_orbital_start(jsite) + 1)%r_cut

          ! Only add the element in case the distances match
          if ( dist <= ir_cut + jr_cut ) call sp%add(io, jo)
        end do
      end do

    end subroutine add_elements

  end subroutine create_sparse_pattern_ac_create

end module esl_create_sparse_pattern_ac_m


