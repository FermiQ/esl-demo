module esl_smear_m
  use prec, only: ip, dp

  implicit none
  private

  public :: smear_t

  !Data structure for the system
  type smear_t

     integer(ip)       :: smearing                     !< which smearing do we use
     real(dp), pointer :: eigenvalues(:,:,:) => null() !< n_state,n_spin,n_kpt
     real(dp)          :: fermi_level
     real(dp)          :: eTemp
     real(dp)          :: eBroad
     real(dp)          :: occ_tol

   contains

     private
     procedure, public :: init
     procedure, public :: calc_fermi_occ
     final :: cleanup

  end type smear_t

  integer, public, parameter :: &
    GAUSSIAN = 0,               &
    FD       = 1,               &
    MP       = 2,               &
    CUBIC    = 3,               &
    COLD     = 4

contains

  subroutine init(this)
    use fdf, only: fdf_get, leqi

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

  subroutine calc_fermi_occ(this, elsic, states)
    use esl_states_m, only: states_t
    use esl_elsi_m, only: elsi_t
 !   use elsi, only: elsi_set_mu_broaden_scheme, &
 !                   elsi_set_mu_broaden_width, &
 !                   elsi_compute_mu_and_occ

    class(smear_t), intent(inout) :: this
    type(elsi_t),   intent(inout) :: elsic
    type(states_t), intent(inout) :: states

  !  call elsi_set_mu_broaden_scheme(elsic%e_h, this%smearing)
  !  call elsi_set_mu_broaden_width(elsic%e_h, this%eBroad)

! TODO: add it back
!    call elsi_compute_mu_and_occ(elsic%e_h, real(states%nel, dp), &
!      & states%nstates, states%nspin, states%nkpt, states%k_weights, &
!      & this%eigenvalues, states%occ_numbers, this%fermi_level)

  end subroutine calc_fermi_occ

end module esl_smear_m
