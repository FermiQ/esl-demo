module esl_mixing_m
  use prec, only : dp,ip

  implicit none
  private

  public :: mixing_t

  ! Data structure for the mixer
  type mixing_t

    real(dp) :: alpha !< Mixing parameter

  contains
    private
    procedure, public :: init
    procedure, public :: linear
    final  :: cleanup
  end type mixing_t

contains

  !Initialize the mixer
  !----------------------------------------------------
  subroutine init(this)
    use fdf, only: fdf_get
    class(mixing_t) :: this

    ! For the moment we read this from SCF.Mix.alpha
    this%alpha = fdf_get('SCF.Mix.alpha', 0.3_dp)

  end subroutine init


  !Release the mixer
  !----------------------------------------------------
  subroutine cleanup(this)
    type(mixing_t) :: this

  end subroutine cleanup

  ! Mix two input vectors
  subroutine linear(this, np, in, out, next)
    class(mixing_t), intent(in) :: this
    integer, intent(in) :: np
    real(dp), intent(in) :: in(:)
    real(dp), intent(in) :: out(:)
    real(dp), intent(out) :: next(:)

    integer :: ip

    do ip = 1, np
      next(ip) = in(ip) * (1._dp - this%alpha) + out(ip) * this%alpha
    end do

  end subroutine linear

end module esl_mixing_m
