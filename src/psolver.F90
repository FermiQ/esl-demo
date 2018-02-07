module psolver_esl
  use Poisson_Solver, only: coulomb_operator, gp
  implicit none
  private

  character(len = 1), parameter :: PSOLVER_3D_PERIOD = 'P'
  character(len = 1), parameter :: PSOLVER_XZ_SURFACE = 'S'
  character(len = 1), parameter :: PSOLVER_Y_WIRE = 'W'
  character(len = 1), parameter :: PSOLVER_ISOLATED = 'F'

  public :: &
       & PSOLVER_3D_PERIOD, &
       & PSOLVER_XZ_SURFACE, &
       & PSOLVER_Y_WIRE, &
       & PSOLVER_ISOLATED

  public :: psolver_t
  
  type psolver_t
     private
     type(coulomb_operator) :: pkernel

     real(gp) :: ionicOffset
     real(gp), dimension(:,:,:), pointer :: ionicPot
   contains
     private
     procedure, public :: init
     procedure, public :: h_potential
     final  :: cleanup
  end type psolver_t

contains

  subroutine init(ps, iproc, nproc, geocode, ndims, hgrids)
    use Poisson_Solver, only: pkernel_init, pkernel_set
    use dictionaries, only: dictionary
    implicit none
    class(psolver_t), intent(out) :: ps
    integer, intent(in) :: iproc, nproc
    character(len = 1), intent(in) :: geocode
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids

    type(dictionary), pointer :: dict_input

    !>@todo Use default options in kernel creation for the moment.
    nullify(dict_input)
    !>@todo Poisson kernel is initialised _and_ computed here, may differ
    !!      the computation later.
    ps%pkernel = pkernel_init(iproc, nproc, dict_input, geocode, ndims, hgrids)
    call pkernel_set(ps%pkernel, verbose=.true.)
    !>@todo Associate the ionic potential later.
    nullify(ps%ionicPot)
    ps%ionicOffset = 0._gp
  end subroutine init

  subroutine cleanup(ps)
    use Poisson_Solver, only: pkernel_free
    implicit none
    type(psolver_t), intent(inout) :: ps

    call pkernel_free(ps%pkernel)
  end subroutine cleanup

  subroutine h_potential(ps, rhopot, ehartree)
    use Poisson_Solver, only: PS_H_potential => H_potential
    implicit none
    class(psolver_t), intent(inout) :: ps
    real(gp), dimension(*), intent(inout) :: rhopot
    real(gp), intent(out) :: ehartree

    call PS_H_potential('G', ps%pkernel, rhopot, ps%ionicPot, ehartree, &
         & ps%ionicOffset, associated(ps%ionicPot))
  end subroutine h_potential
end module psolver_esl
