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
    use esl_basis_m, only: PLANEWAVES, ATOMCENTERED
    
    type(system_t), intent(inout) :: system

    ! Move quantities required for a move

    select case ( system%basis%type )
    case ( PLANEWAVES )
      call next_planewave()
    case ( ATOMCENTERED )
      call next_atomicorbs()
    end select

  contains

    subroutine next_atomicorbs()

      ! add content for the initialization for a new plane-wave step
      call system%update(periodic=.false.)

    end subroutine next_atomicorbs

    subroutine next_planewave()

      ! add content for the initialization for a new plane-wave step
      call system%update(periodic=.false.)

    end subroutine next_planewave

  end subroutine next_step_setup

end module esl_next_step_m

