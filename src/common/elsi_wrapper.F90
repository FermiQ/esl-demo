!> This module interfaces with ELSI.
!> TODO: Read input options
module esl_elsi_m

  use prec, only: ip, dp
  use elsi, only: elsi_handle, elsi_init, elsi_finalize 

  implicit none
  private

  public :: elsi_t

  !< Data structure for ELSI.
  type elsi_t
    type(elsi_handle) :: e_h
    real(dp)          :: KS_energy
    real(dp)          :: fermi_level
    real(dp)          :: entropy
  contains
    private
    procedure, public :: init
    final :: cleanup
  end type elsi_t

  integer(ip), parameter :: ELPA       = 1
  integer(ip), parameter :: MULTI_PROC = 1
  integer(ip), parameter :: SIESTA_CSR = 2

contains

  !> Initialize ELSI
  subroutine init(this, n_basis, n_electron, n_state)

    class(elsi_t), intent(inout) :: this
    integer(ip),   intent(in)    :: n_basis
    real(dp),      intent(in)    :: n_electron
    integer(ip),   intent(in)    :: n_state

    ! Initialize an ELSI handle
    call elsi_init(this%e_h, ELPA, MULTI_PROC, SIESTA_CSR, n_basis, &
      & n_electron, min(2*int(n_electron), n_basis))

    this%KS_energy   = 0._dp
    this%fermi_level = 0._dp
    this%entropy     = 0._dp

  end subroutine init

  !> Finalize ELSI
  subroutine cleanup(this)
    type(elsi_t), intent(inout) :: this

    call elsi_finalize(this%e_h)

  end subroutine cleanup

end module esl_elsi_m
