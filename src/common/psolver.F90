module esl_psolver_m
  use dictionaries
  use Poisson_Solver
  implicit none
  private

  character(len=1), parameter :: PSOLVER_3D_PERIOD = 'P'
  character(len=1), parameter :: PSOLVER_XZ_SURFACE = 'S'
  character(len=1), parameter :: PSOLVER_Y_WIRE = 'W'
  character(len=1), parameter :: PSOLVER_ISOLATED = 'F'

  public :: &
      & PSOLVER_3D_PERIOD, &
      & PSOLVER_XZ_SURFACE, &
      & PSOLVER_Y_WIRE, &
      & PSOLVER_ISOLATED

  public :: psolver_t

  type psolver_t
    private
    type(coulomb_operator) :: pkernel
  contains
    procedure, public :: init
    procedure, public :: hartree_potential
    final  :: finalizer
  end type psolver_t

contains

  subroutine init(ps, iproc, nproc, geocode, ndims, hgrids)
    class(psolver_t),   intent(inout) :: ps
    integer,            intent(in)    :: iproc, nproc
    character(len = 1), intent(in)    :: geocode
    integer,            intent(in)    :: ndims(3)
    real(gp),           intent(in)    :: hgrids(3)

    type(dictionary), pointer :: dict_input

    !>@todo Use default options in kernel creation for the moment.
    nullify(dict_input)
    call dict_init(dict_input)
    !>@todo Poisson kernel is initialised _and_ computed here, may differ
    !!      the computation later.
    write(*,*) ndims, hgrids
    ps%pkernel = pkernel_init(iproc, nproc, dict_input, geocode, ndims, hgrids)
    call dict_free(dict_input)
    call pkernel_set(ps%pkernel, verbose=.true.)
  end subroutine init

  subroutine finalizer(ps)
    type(psolver_t), intent(inout) :: ps

    call pkernel_free(ps%pkernel)
  end subroutine finalizer

  subroutine hartree_potential(ps, hartree, ionicPot, ionicOffset, ehartree)
    class(psolver_t), intent(inout) :: ps
    real(gp),         intent(inout) :: hartree(*) !< On input should hold the density. On output it will contain the Hartree potential
    real(gp),         intent(inout) :: ionicPot(*)
    real(gp),         intent(in)    :: ionicOffset
    real(gp),         intent(out)   :: ehartree

    call H_potential('G', ps%pkernel, hartree, ionicPot, ehartree, ionicOffset, .true., quiet = "YES")

  end subroutine hartree_potential

end module esl_psolver_m
