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
  subroutine create_sparse_pattern_ac_create(sys, sp)
    
    use prec, only: dp
    use system_esl, only: system_t
    use sparse_pattern, only: sparse_pattern_t
    
    class(system_t), intent(in) :: sys
    type(sparse_pattern_t), intent(inout) :: sp
    integer :: no, max_no, io, jo
    integer :: ia, ja
    real(dp) :: r2, dist

    ! Start by deallocation of the sparse pattern
    call sp%delete()

    ! Figure out the number of orbitals
    no = 0
    max_no = 0
    do ia = 1, sys%nAtoms
       max_no = max(max_no, 1) ! 1 should be replaced by number of orbitals per atom
       no = no + 1 ! to be filled
    end do

    ! Now re-initialize the sparse matrix.
    ! In this case we will assume a maximum of 20 atomic connections
    call sp%init(no, no, np=max_no * 20)

    ! Loop over all atoms
    do ia = 1, sys%nAtoms - 1
       is = sys%ispecie(ia)

       ! Add the connections to it-self
       call add_elements(ia, ia, 0._dp)

       ! Only loop the remaining atoms (no need to double process)
       do ja = ia + 1, sys%nAtoms
          js = sys%ispecie(ja)

          ! Calculate whether the distance between the two
          ! atoms is within their basis range.
          r2 = sys%pseudo(is)%rmax + sys%pseudo(js)%rmax

          ! Calculate the distance between the two atomic centers
          dist = sqrt(sum((sys%xyz(:,ia) - sys%xyz(:,ja)) ** 2))

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
      do io = sys%first_orb(ia) , sys%first_orb(ia + 1) - 1
         do jo = sys%first_orb(ja) , sys%first_orb(ja + 1) - 1
            call sp%add(io, jo)
         end do
      end do
      
    end subroutine add_elements

  end subroutine create_sparse_pattern_ac_create
  
end module esl_create_sparse_pattern_ac_m
    

