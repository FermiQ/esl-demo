module esl_species_m

  use prec, only: dp
  use esl_message_m
  use pspiof_m

  implicit none
  private

  public :: species_t

  !Data structure for the pseudos
  type species_t
    character(len=10) :: label
    type(pspiof_pspdata_t) :: psp

    ! TODO fix this argument
    real(dp) :: Z = 0._dp !< Ionic charge, real because of mixed species
  contains
    private
    procedure, public :: init
    final :: cleanup
  end type species_t

contains

  !Initialize the species data
  !----------------------------------------------------
  subroutine init(this, label, filename)
    class(species_t) :: this
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: filename

    integer :: ierr

    this%label = label
    
    ! Make sure pseudopotential data can be allocated
    ierr = check_error_pspio(pspiof_pspdata_alloc(this%psp))
    if ( ierr /= PSPIO_SUCCESS ) return

    ! Parse the file and record the format
    !ierr = check_error_pspio(pspiof_pspdata_read(this%psp, PSPIO_FMT_UNKNOWN, filename))
    if ( ierr /= PSPIO_SUCCESS ) then
       call message_error("Could not read pseudopotential from file: "//filename)
       return
    endif

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(species_t) :: this

    call pspiof_pspdata_free(this%psp)

  end subroutine cleanup

  !----------------------------------------------------
  !Private routines
  !----------------------------------------------------

  !Check error codes from Pspio
  !----------------------------------------------------
  function check_error_pspio(retcode) result(ierr)
    integer, intent(in) :: retcode

    integer :: ierr

    ierr = retcode
    if ( retcode /= PSPIO_SUCCESS ) then
       call pspiof_error_flush()
    end if

  end function check_error_pspio

end module esl_species_m
