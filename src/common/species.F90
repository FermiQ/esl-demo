module esl_species_m
  use esl_grid_m
  use esl_message_m
  use prec
  use pspiof_m
  use yaml_output
  implicit none
  private

  public :: species_t

  !Data structure for the pseudos
  type species_t
    character(len=10) :: label
    
    type(pspiof_pspdata_t)  :: psp
    type(pspiof_meshfunc_t) :: rho

    integer :: n_orbitals
    type(pspiof_state_t), allocatable :: orbitals(:)

    real(dp) :: z_ion = 0._dp !< Ionic charge, real because of mixed species
    real(dp) :: q = 0._dp     !< Electronic charge
  contains
    private
    procedure, public :: init
    procedure, public :: get_orbital
    procedure, public :: summary
    final :: cleanup
  end type species_t

contains

  !Initialize the species data
  !----------------------------------------------------
  subroutine init(this, label, filename)
    class(species_t) :: this
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: filename

    integer :: ierr, io

    this%label = label
    
    ! Make sure pseudopotential data can be allocated
    ierr = check_error_pspio(pspiof_pspdata_alloc(this%psp))
    if ( ierr /= PSPIO_SUCCESS ) return

    ! Parse the file and record the format
    ierr = check_error_pspio(pspiof_pspdata_read(this%psp, PSPIO_FMT_UNKNOWN, trim(filename)))
    if ( ierr /= PSPIO_SUCCESS ) then
       call message_error("Could not read pseudopotential from file: "//filename)
       return
    endif

    ! Store some information
    this%rho = pspiof_pspdata_get_rho_valence(this%psp)
    this%n_orbitals = pspiof_pspdata_get_n_states(this%psp)
    if (this%n_orbitals > 0) then
      allocate(this%orbitals(this%n_orbitals))
      do io = 1, this%n_orbitals
        this%orbitals(io) = pspiof_pspdata_get_state(this%psp, io)
      end do
    end if
    
    this%z_ion = pspiof_pspdata_get_zvalence(this%psp)
    this%q = pspiof_pspdata_get_nelvalence(this%psp)

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(species_t) :: this

    integer :: io
    
    call pspiof_pspdata_free(this%psp)
    call pspiof_meshfunc_free(this%rho)
    if (allocated(this%orbitals)) then
      do io = 1, this%n_orbitals
        call pspiof_state_free(this%orbitals(io))
      end do
    end if
    
  end subroutine cleanup

  !summary
  !----------------------------------------------------
  subroutine summary(this)
    class(species_t) :: this

    call yaml_mapping_open(this%label)
    call yaml_map("Ionic charge", this%z_ion)
    call yaml_map("Electronic charge", this%q)
    call yaml_mapping_close()

  end subroutine summary

  
  subroutine get_orbital(this, io, orbital, ll, occ)
    class(species_t) :: this
    integer,                 intent(in)  :: io
    type(pspiof_meshfunc_t), intent(out) :: orbital
    integer,                 intent(out) :: ll
    real(dp),                intent(out) :: occ
    
    orbital = pspiof_state_get_wf(this%orbitals(io))
    ll = pspiof_qn_get_l(pspiof_state_get_qn(this%orbitals(io)))
    occ = pspiof_state_get_occ(this%orbitals(io))

  end subroutine get_orbital
  
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
!       call pspiof_error_flush()
    end if

  end function check_error_pspio

end module esl_species_m
