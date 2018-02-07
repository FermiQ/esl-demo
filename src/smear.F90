module smear_esl
  use prec, only : dp,ip
  use elsi_wrapper_esl, only: elsi_calc_fermi_and_occ
  use fdf, only  : fdf_double, fdf_string

  implicit none
  private

  public :: smear_t
  public :: smear_calc_fermi_and_occ

  !Data structure for the system
  type smear_t
    integer(kind=ip)           :: smearing           !< which smearing do we use
    real(kind=dp), pointer     :: eigenvalues(:,:,:) !< n_state,n_spin,n_kpt
    real(kind=dp), allocatable :: occ_numbers(:,:,:)
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
    end cleanup

    !Compute the Fermi level and occupation numbers
    !----------------------------------------------------
    subroutine smear_calc_fermi_and_occ(this)
      type(smear_t), intent(inout) :: this

      call elsi_calc_fermi_and_occ(this%smearing,this%eigenvalues,&
        this%occ_numbers,this%fermi_level)

    end subroutine smear_calc_fermi_and_occ

  end module smear_esl
