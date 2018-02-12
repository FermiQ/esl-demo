module esl_smear_m
  use prec, only : ip, dp
  use elsi_wrapper_esl, only: elsi_calc_fermi_and_occ
  use fdf, only  : fdf_get, leqi

  use esl_states_m

  implicit none
  private

  public :: smear_t
  public :: smear_calc_fermi_and_occ

  !Data structure for the system
  type smear_t

     integer(ip)           :: smearing           !< which smearing do we use
     real(dp), pointer     :: eigenvalues(:,:,:) => null() !< n_state,n_spin,n_kpt
     real(dp)              :: fermi_level
     real(dp)              :: eTemp
     real(dp)              :: eBroad
     real(dp)              :: occ_tol

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

    sm = fdf_get('Smearing','Gaussian')

    if ( leqi(sm, 'Gaussian') )  Then
       this%smearing = GAUSSIAN
    else if ( leqi(sm, 'FD') ) then
       this%smearing = FD
    else if ( leqi(sm, 'MP') ) then
       this%smearing = MP
    else if ( leqi(sm, 'CUBIC') ) then
       this%smearing = CUBIC
    else if (leqi(sm, 'COLD') ) then
       this%smearing = COLD
    else
    endif

    ! 0.00095 Ha ~ 300 K
    this%eTemp = fdf_get('Smearing.Temp', 0.00095_dp, 'Ha')
    this%eBroad = fdf_get('eBroad', 0.1_dp, 'Ha')

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

    if ( this%smearing /= COLD ) then
       call elsi_calc_fermi_and_occ(elsi,states%nel, states%nstates, states%nspin, states%nkpt, &
            & this%eigenvalues, states%occ_numbers, states%k_weights, this%fermi_level)
    end if

  end subroutine smear_calc_fermi_and_occ

end module esl_smear_m
