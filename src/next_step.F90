!< Subroutine used to start a new step around an SCF loop
!< This could be one of (but not limited to):
!<  1. MD step
!<  2. change of cutoff values
!<  3. change of XXX
module esl_next_step_m

  implicit none

  public :: next_step_setup

contains

  subroutine next_step_setup(system)

    use esl_system_m, only: system_t
    use esl_basis_m, only: PLANEWAVES, ATOMICORBS

    !< System that we wish to process as a new step
    class(system_t), intent(inout) :: system

    ! TODO insert generic initialization stuff

    select case ( system%basis%type )
    case ( PLANEWAVES )
       call next_planewave()
    case ( ATOMICORBS )
       call next_atomicorbs()
    end select

  contains

    subroutine next_atomicorbs()

      use esl_sparse_pattern_m, only: sparse_pattern_t
      use esl_create_sparse_pattern_ac_m, only: create_sparse_pattern_ac_create
      use esl_overlap_matrix_ac_m, only: overlap_matrix_ac_calculate
      use esl_density_matrix_ac_m, only: density_matrix_ac_next

      ! Target is necessary to preserve old pointers
      type(sparse_pattern_t), target :: old_sp

      ! Preserve the old sparse pattern such that the
      ! sparse matrices still points to the pattern.
      call move_alloc(system%sparse_pattern_t, old_sp)

      ! copy old sparse pattern
      ! Initialize the sparse pattern
      call create_sparse_pattern_ac_create(system, system%sparse_pattern)

      ! Essentially we have to figure out whether the previous
      ! sparse patterns and quantities needs mangling (i.e.
      ! copy the "old" DM -> "new" DM.
      ! There are various choices here.
      ! TODO Handle new-step DM changes

      ! Calculate the overlap matrix
      ! This requires that the grid information is present
      call overlap_matrix_ac_calculate(system, &
           system%sparse_pattern, system%S)

      ! Figure out the next density matrix
      ! This routine will initially fill it with the
      ! atomic fillings.
      call density_matrix_ac_next(system, old_sp, &
           system%sparse_pattern, system%DM)

      ! Clean-up the old sparse-matrix
      call old_sp%delete()

    end subroutine next_atomicorbs

    subroutine next_planewave()

      ! add content for the initialization for a new plane-wave step

    end subroutine next_planewave

  end subroutine next_step_setup

end module esl_next_step_m

