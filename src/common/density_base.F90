module esl_density_base_m
  use prec, only : dp
  
  implicit none

  private

  public :: density_base_t

  type density_base_t
    private
    real(dp), allocatable, public :: rho(:)
  end type density_base_t

end module esl_density_base_m
