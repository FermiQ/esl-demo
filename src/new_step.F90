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

      ! Initialize the sparse pattern
      call init_sparse_pattern(system, system%sparse_pattern)

    end subroutine new_atomicorbs

    subroutine new_planewave()

      ! add content for the initialization for a new plane-wave step
      
    end subroutine new_planewave

  end subroutine new_step

end module new_step_esl
  
