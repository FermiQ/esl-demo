module smear_esl
  use elsi_wrapper_esl, only :: elsi_calc_fermi_and_occ

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
  end type smear_t

   integer, public, parameter :: &
    GAUSSIAN = 0,           &
    FD       = 1,           &
    MP       = 3,           &
    CUBIC    = 4,           &
    COLD     = 5

  contains

   !Compute the Fermi level and occupation numbers
   !----------------------------------------------------
   subroutine smear_calc_fermi_and_occ(this)
     type(smear_t), intent(inout) :: this

     call elsi_calc_fermi_and_occ(this%smearing,this%eigenvalues,&
          this%occ_numbers,this%fermi_level)

   end subroutine smear_calc_fermi_and_occ

end module smear_esl
