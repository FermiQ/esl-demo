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

    use esl_basis_ac_m, only: basis_ac_t
    use esl_sparse_pattern_m, only: sparse_pattern_t
    use esl_sparse_matrix_m, only: sparse_matrix_t

    class(basis_ac_t), intent(inout) :: basis
    class(sparse_pattern_t), intent(inout) :: sp
    class(sparse_matrix_t), intent(inout) :: S

    integer :: ia, is, io, iio, ind, jo, ja, js, jjo

    ! Re-initialize the sparse matrix
    call S%init(sp)

    ! Loop over all orbital connections in the sparse pattern and
    ! calculate the overlap matrix for each of them
    do ia = 1, basis%n_sites
      is = basis%species_idx(ia)

      ! Loop on orbitals
      do io = basis%site_function_start(ia), basis%site_function_start(ia + 1) - 1
        ! Orbital index on atom
        iio = io - basis%site_function_start(ia) + 1

        ! Loop entries in the sparse pattern
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

          ! Figure out which atom this orbital belongs too
          jo = sp%column(ind)
          ! Figure out the atomic index of the orbital
          ja = basis%function_site(jo)
          js = basis%species_idx(ja)
          jjo = jo - basis%site_function_start(ja) + 1

          ! We are now in a position to calculate the
          ! overlap matrix. I.e. we know the atom, the
          ! orbital indices and their positions

          ! TODO fix calls here
!          S%M(ind) = &
!              basis%grid%overlap(xyz(:,ia), TODO-AO, pseudo(is)%rmax, &
!              xyz(:,ja), TODO-AO, pseudo(js)%rmax)

        end do

      end do

    end do

  end subroutine overlap_matrix_ac_calculate

end module esl_overlap_matrix_ac_m
