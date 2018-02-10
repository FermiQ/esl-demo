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
  subroutine create_sparse_pattern_ac_create(geo, basis, sp)

    use prec, only: dp
    use esl_geometry_m, only: geometry_t
    use esl_basis_ac_m, only: basis_ac_t
    use esl_sparse_pattern_m, only: sparse_pattern_t

    class(geometry_t), intent(in) :: geo
    class(basis_ac_t), intent(in) :: basis
    type(sparse_pattern_t), intent(inout) :: sp
    integer :: no, max_no, io, jo
    integer :: ia, ja, is, js
    real(dp) :: r2, dist

    ! Start by deallocation of the sparse pattern
    call sp%delete()

    ! Figure out the number of orbitals
    no = 0
    max_no = 0
    do ia = 1, geo%n_atoms
      max_no = max(max_no, 1) ! 1 should be replaced by number of orbitals per atom
      no = no + 1 ! to be filled
    end do

    ! Now re-initialize the sparse matrix.
    ! In this case we will assume a maximum of 20 atomic connections
    call sp%init(no, no, np=max_no * 20)

    ! Loop over all atoms
    do ia = 1, geo%n_atoms
      is = basis%species_idx(ia)

      ! Add the connections to it-self
      call add_elements(ia, ia, 0._dp)

      ! Only loop the remaining atoms (no need to double process)
      do ja = ia + 1, geo%n_atoms
        js = basis%species_idx(ja)

        ! Calculate whether the distance between the two
        ! atoms is within their basis range.
! TODO FIX RMAX         
!          r2 = pseudo(is)%rmax + pseudo(js)%rmax

        ! Calculate the distance between the two atomic centers
        dist = sqrt(sum((geo%xyz(:,ia) - geo%xyz(:,ja)) ** 2))

        ! Only process if the maximum distance is within range.
        if ( dist <= r2 ) then

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
      integer :: io, jo

      ! TODO do orbital dependent distances

      ! Loop orbitals on both atoms
      do io = basis%site_function_start(ia) , basis%site_function_start(ia+1) - 1
        do jo = basis%site_function_start(ja) , basis%site_function_start(ja+1) - 1
          call sp%add(io, jo)
        end do
      end do

    end subroutine add_elements

  end subroutine create_sparse_pattern_ac_create

end module esl_create_sparse_pattern_ac_m


