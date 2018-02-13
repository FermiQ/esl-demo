!< This routine handles the calculation of the overlap matrix
!<
!< It expects a pre-allocated sparse pattern which is used as the
!< pre-cursor for when two orbitals are overlapping.
!< This also allows one to "play" with shrinking sparse matrices
!< in certain cases.
module esl_overlap_matrix_ac_m

  implicit none

  public :: overlap_matrix_ac_calculate

contains

  subroutine overlap_matrix_ac_calculate(basis, sp, S)
    use prec, only: dp
    use esl_basis_ac_m, only: basis_ac_t
    use esl_sparse_pattern_m, only: sparse_pattern_t
    use esl_sparse_matrix_m, only: sparse_matrix_t

    class(basis_ac_t), intent(inout) :: basis
    class(sparse_pattern_t), intent(inout) :: sp
    class(sparse_matrix_t), intent(inout) :: S

    integer :: ia, is, io, iio, ind, jo, ja, js, jjo
    real(kind=dp), allocatable :: ao1(:), ao2(:)
    integer :: ll, mm

!    allocate(ao1(1:basis%grid%np))
!    allocate(ao2(1:basis%grid%np))

    ! Re-initialize the sparse matrix
    call S%init(sp)

    ! Loop over all orbital connections in the sparse pattern and
    ! calculate the overlap matrix for each of them
    do ia = 1, basis%n_site
      is = basis%site_state_idx(ia)

      ! Loop on orbitals
      do io = basis%site_orbital_start(ia), basis%site_orbital_start(ia + 1) - 1
        ! Orbital index on atom
        iio = io - basis%site_orbital_start(ia) + 1

        !TODO: get ll and mm
 !       call basis%grid%radial_function(basis%orbitals(iio), ll, mm, xyz(:,ia), ao1)

        ! Loop entries in the sparse pattern
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

          ! Figure out which atom this orbital belongs too
          jo = sp%column(ind)
          ! Figure out the atomic index of the orbital
          ja = basis%orbital_site(jo)
          js = basis%site_state_idx(ja)
          jjo = jo - basis%site_orbital_start(ja) + 1

          ! We are now in a position to calculate the
          ! overlap matrix. I.e. we know the atom, the
          ! orbital indices and their positions

          !TODO: get ll and mm
!          call basis%grid%radial_function(basis%orbitals(iio), ll, mm, xyz(:,ja), ao2)

!          S%M(ind) = &
!              basis%grid%overlap(xyz(:,ia), ao1, pseudo(is)%rmax, &
!              xyz(:,ja), ao2, pseudo(js)%rmax)

        end do

      end do

    end do

    deallocate(ao1) 
    deallocate(ao2)

  end subroutine overlap_matrix_ac_calculate

end module esl_overlap_matrix_ac_m
