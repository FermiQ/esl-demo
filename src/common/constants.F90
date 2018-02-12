module esl_constants_m

  use prec, only: dp

  implicit none

  ! privatize all imported variables
  private :: dp

  public

  real(dp), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter :: DEG = PI / 180._dp

end module esl_constants_m
