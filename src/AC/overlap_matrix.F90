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

  subroutine overlap_matrix_ac_calculate(basis, grid, sp, S)
    use prec, only: dp
    use esl_basis_ac_m, only: basis_ac_t
    use esl_grid_m, only: grid_t
    use esl_sparse_pattern_m, only: sparse_pattern_t
    use esl_sparse_matrix_m, only: sparse_matrix_t

    class(basis_ac_t), intent(inout) :: basis
    type(grid_t), intent(in) :: grid
    type(sparse_pattern_t), intent(in) :: sp
    type(sparse_matrix_t), intent(inout) :: S

    integer :: ia, is, io, iio, ind, jo, ja, js, jjo
    real(dp), allocatable :: iao(:), jao(:)
    real(dp) :: ixyz(3), ir_max, jxyz(3), jr_max
    integer :: il, im, jl, jm

    allocate(iao(1:grid%np))
    allocate(jao(1:grid%np))

    ! Re-initialize the sparse matrix
    call S%init(sp)

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
        call grid%radial_function_ylm(basis%state(is)%orb(iio)%R, il, im, ixyz(:), iao)

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
          call grid%radial_function_ylm(basis%state(js)%orb(jjo)%R, jl, jm, jxyz(:), jao)

          S%M(ind) = &
              grid%overlap(ixyz(:), iao, ir_max, jxyz(:), jao, jr_max)

          ! DEBUG print
          if ( ia == ja .and. iio == jjo ) &
              print *,' Diagonal overlap matrix: ', ia, iio, S%M(ind)
          !print *,' Calculing overlap matrix: ', ia, iio, ja, jjo, S%M(ind)

        end do

      end do

    end do

    deallocate(iao,jao)

  end subroutine overlap_matrix_ac_calculate

end module esl_overlap_matrix_ac_m
