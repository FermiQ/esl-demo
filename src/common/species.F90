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
    type(pspiof_meshfunc_t) :: vlocal

    integer :: n_orbitals
    type(pspiof_state_t), allocatable :: orbitals(:)

    real(dp) :: z_ion = 0._dp !< Ionic charge, real because of mixed species
    real(dp) :: q = 0._dp     !< Electronic charge
  contains
    private
    procedure, public :: init
    procedure, public :: get_radial_orbital
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

    integer :: ierr, io, ir
    real(dp), pointer :: r(:)
    real(dp), allocatable :: vl(:)
    type(pspiof_potential_t) :: vlocal
    type(pspiof_mesh_t) :: mesh

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

    ! Get local potential.
    ! At the moment pspio does not have a getter for the meshfunc_t from a potential_t, so we are going to use a workaround
    mesh = pspiof_pspdata_get_mesh(this%psp)
    vlocal = pspiof_pspdata_get_vlocal(this%psp)
    allocate(vl(pspiof_mesh_get_np(mesh)))
    ierr = check_error_pspio(pspiof_meshfunc_alloc(this%vlocal, pspiof_mesh_get_np(mesh)))
    r => pspiof_mesh_get_r(mesh)
    do ir = 1, pspiof_mesh_get_np(mesh)
      vl(ir) = pspiof_potential_eval(vlocal, r(ir))
    end do
    ierr = check_error_pspio(pspiof_meshfunc_init(this%vlocal, mesh, vl))
    deallocate(vl)

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
    call pspiof_meshfunc_free(this%vlocal)

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

  
  subroutine get_radial_orbital(this, io, ll, orbital, occ)
    class(species_t) :: this
    integer,                           intent(in)  :: io
    integer,                 optional, intent(out) :: ll
    type(pspiof_meshfunc_t), optional, intent(out) :: orbital
    real(dp),                optional, intent(out) :: occ

    if (present(ll)) then
      ll = pspiof_qn_get_l(pspiof_state_get_qn(this%orbitals(io)))
    end if
    if (present(orbital)) then
      orbital = pspiof_state_get_wf(this%orbitals(io))
    end if
    if (present(occ)) then
      occ = pspiof_state_get_occ(this%orbitals(io))
    end if

  end subroutine get_radial_orbital
  
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
