module esl_basis_base_m
  use esl_grid_m
  
  implicit none

  private

  public :: basis_base_t

  type basis_base_t
    private
    integer,      public :: size
    type(grid_t), public :: grid !< Auxiliary real-space grid
  end type basis_base_t

end module esl_basis_base_m
