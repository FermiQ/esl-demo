!< Module for handling density matrices
!<
!< This module should implement the basics of creating
!< an initial (valence filled) density matrix and/or
!< read-in/extrapolate a DM to a different structure.
module esl_density_matrix_ac_m

  use esl_system_m, only: system_t
  use esl_sparse_pattern_m, only: sparse_pattern_t
  use esl_sparse_matrix_m, only: sparse_matrix_t

  implicit none

  public :: density_matrix_ac_next
  public :: density_matrix_ac_init_atomic

contains

  !< Prepare the next density matrix
  !<
  !< This routine gets passed an old sparse pattern and a
  !< new sparse pattern.
  !< When the old sparse pattern is not allocated (created)
  !< we automatically initialize the DM with the atomic fillings.
  subroutine density_matrix_ac_next(sys, old_sp, new_sp, DM)

    class(system_t), intent(in) :: sys
    type(sparse_pattern_t), intent(in) :: old_sp
    type(sparse_pattern_t), intent(in), target :: new_sp
    ! One element per spin component.
    type(sparse_matrix_t), intent(inout), allocatable :: DM(:)

    if ( old_sp%initialized() ) then

      ! For now we still do the atomic fillings...
      call density_matrix_ac_init_atomic(sys, new_sp, DM)

    else

      call density_matrix_ac_init_atomic(sys, new_sp, DM)

    end if

  end subroutine density_matrix_ac_next

  !< Initialize the diagonal density matrix with atomic fillings
  subroutine density_matrix_ac_init_atomic(sys, sp, DM)

    use prec, only: dp

    class(system_t), intent(in) :: sys
    type(sparse_pattern_t), intent(in), target :: sp
    ! One element per spin component.
    type(sparse_matrix_t), intent(inout), allocatable :: DM(:)

    real(dp) :: frac_s
    integer :: nspin
    integer :: ia, is, io, iio, ind

    nspin = 1
    ! Initialize the density matrix
    if ( allocated(DM) ) then
      ! We assume it is allocated correctly (nspin is correct)
    else
      ! TODO insert correct size  (number of spin-components)
      allocate(DM(nspin))
    end if

    do is = 1, nspin
      call DM(is)%init(sp)

      ! Initialize all elements to zero
      DM(is)%M = 0._dp
    end do

    ! Define the fraction of filling per spin channel
    frac_s = 1._dp / nspin

    ! Loop over all orbital connections in the sparse pattern and
    ! set the diagonal density matrix
    do ia = 1, sys%geo%n_atoms
      is = sys%geo%species_idx(ia)

      ! Loop on orbitals
      do io = sys%basis%ac%site_function_start(ia), sys%basis%ac%site_function_start(ia + 1) - 1
        ! Orbital index on atom
        iio = io - sys%basis%ac%site_function_start(ia) + 1

        ! Loop entries in the sparse pattern
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1
          ! Skip if it does not correspond to the diagonal element
          if ( sp%column(ind) /= io ) cycle

          ! Set diagonal component
! TODO add atomic charge specification
!          DM(:)%M(ind) = sys%basis% TODO basis %info(is)%q0(iio) * frac_s

        end do

      end do

    end do

  end subroutine density_matrix_ac_init_atomic

end module esl_density_matrix_ac_m
