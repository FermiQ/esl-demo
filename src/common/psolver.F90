module esl_psolver_m

  use Poisson_Solver, only: coulomb_operator, gp
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
    procedure, public :: h_potential
    final  :: cleanup
  end type psolver_t

contains

  subroutine init(ps, iproc, nproc, geocode, ndims, hgrids)
    use Poisson_Solver, only: pkernel_init, pkernel_set
    use dictionaries, only: dictionary, dict_init, dict_free

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
    ps%pkernel = pkernel_init(iproc, nproc, dict_input, geocode, ndims, hgrids)
    call dict_free(dict_input)
    call pkernel_set(ps%pkernel, verbose=.true.)
  end subroutine init

  subroutine cleanup(ps)
    use Poisson_Solver, only: pkernel_free
    type(psolver_t), intent(inout) :: ps

    call pkernel_free(ps%pkernel)
  end subroutine cleanup

  subroutine h_potential(ps, hartree, ionicPot, ionicOffset, ehartree)
    use Poisson_Solver, only: PS_H_potential => H_potential

    class(psolver_t), intent(inout) :: ps
    real(gp),         intent(inout) :: hartree(*) !< On input should hold the density. On output it will contain the Hartree potential
    real(gp),         intent(inout) :: ionicPot(*)
    real(gp),         intent(in)    :: ionicOffset
    real(gp),         intent(out)   :: ehartree

    call PS_H_potential('G', ps%pkernel, hartree, ionicPot, ehartree, ionicOffset, .true., quiet = "YES")

  end subroutine h_potential

end module esl_psolver_m
