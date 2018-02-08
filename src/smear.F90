module smear_esl
  use prec, only : ip, dp
  use elsi_wrapper_esl, only: elsi_calc_fermi_and_occ
  use fdf, only  : fdf_double, fdf_string

  use states_esl

  implicit none
  private

  public :: smear_t
  public :: smear_calc_fermi_and_occ

  !Data structure for the system
  type smear_t
    integer(kind=ip)           :: smearing           !< which smearing do we use
    real(kind=dp), pointer     :: eigenvalues(:,:,:) !< n_state,n_spin,n_kpt
    real(kind=dp)              :: fermi_level
    real(kind=dp)              :: eTemp
    real(kind=dp)              :: eBroad
    real(kind=dp)              :: occ_tol
  contains
    private
    procedure, public :: init
    final :: cleanup
  end type smear_t

  integer, public, parameter :: &
    GAUSSIAN = 0,           &
    FD       = 1,           &
    MP       = 3,           &
    CUBIC    = 4,           &
    COLD     = 5

contains

  subroutine init(this)
  class(smear_t) :: this

    character(len=100) :: sm
    sm=fdf_string('Smearing','Gaussian')
    if (sm == 'Gaussian')  Then
      this%smearing=GAUSSIAN
    else if (sm=='FD') then
      this%smearing=FD
    else if (sm == 'MP') then
      this%smearing=MP
    else if (sm == 'CUBIC') then
      this%smearing=CUBIC
    else if (sm=='COLD') then
      this%smearing = COLD
    else

    endif 
    this%eTemp=fdf_double('eTemp', 300.0_dp)
    this%eBroad=fdf_double('eBroad',0.1_dp)
  end subroutine init

  subroutine cleanup(this)
    type(smear_t) :: this
  end subroutine cleanup

   !Compute the Fermi level and occupation numbers
   !----------------------------------------------------
  subroutine smear_calc_fermi_and_occ(this, elsi, states)
    use elsi_wrapper_esl
     type(smear_t), intent(inout) :: this
     type(elsi_t), intent(inout) :: elsi
     type(states_t),   intent(in) :: states

     integer :: n_electron
     !@todo: put this in states
     real(dp), dimension(states%nkpt) :: k_weights

     !@todo: don't compute n_electron
     n_electron = sum(states%occ_numbers(:,:,1))
     call elsi_calc_fermi_and_occ(elsi, n_electron, states%nstates, states%nspin, states%nkpt, &
          & this%eigenvalues, states%occ_numbers, k_weights, this%fermi_level)

    end subroutine smear_calc_fermi_and_occ

  end module smear_esl
