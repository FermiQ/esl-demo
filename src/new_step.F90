!< Subroutine used to start a new step around an SCF loop
!< This could be one of (but not limited to):
!<  1. MD step
!<  2. change of cutoff values
!<  3. change of XXX
module new_step_esl

  implicit none

  public :: new_step

contains

  subroutine new_step(system)

    use system_esl, only: system_t
    use basis_esl, only: PLANEWAVES, ATOMICORBS

    !< System that we wish to process as a new step
    class(system_t), intent(inout) :: system

    ! TODO insert generic initialization stuff

    select case ( system%basis%type )
    case ( PLANEWAVES )
       call new_planewave()
    case ( ATOMICORBS )
       call new_atomicorbs()
    end select

  contains

    subroutine new_atomicorbs()
      
      use init_sparse_pattern_esl, only: init_sparse_pattern
      use overlap_matrix_esl, only: calc_overlap_matrix
      use density_matrix_esl, only: next_density_matrix

      ! Target is necessary to preserve old pointers
      type(sparse_pattern_t), target :: old_sp

      ! Preserve the old sparse pattern such that the
      ! sparse matrices still points to the pattern.
      call move_alloc(system%sparse_pattern_t, old_sp)
      
      ! copy old sparse pattern
      ! Initialize the sparse pattern
      call init_sparse_pattern(system, system%sparse_pattern)

      ! Essentially we have to figure out whether the previous
      ! sparse patterns and quantities needs mangling (i.e.
      ! copy the "old" DM -> "new" DM.
      ! There are various choices here.
      ! TODO Handle new-step DM changes

      ! Calculate the overlap matrix
      ! This requires that the grid information is present
      call calc_overlap_matrix(system, &
           system%sparse_pattern, system%S)

      ! Figure out the next density matrix
      ! This routine will initially fill it with the
      ! atomic fillings.
      call next_density_matrix(system, old_sp, &
           system%sparse_pattern, system%DM)

      ! Clean-up the old sparse-matrix
      call old_sp%delete()

    end subroutine new_atomicorbs

    subroutine new_planewave()

      ! add content for the initialization for a new plane-wave step
      
    end subroutine new_planewave

  end subroutine new_step

end module new_step_esl
  
