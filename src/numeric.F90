module numeric_esl
  use prec, only : dp,ip
  use fdf, only : fdf_integer

  implicit none
  private

  public ::                 &
            invert_cell, &
            init_random

contains

  !----------------------------------------------------
  pure real(kind=dp) function determinant(a)
    real(kind=dp), intent(in)  :: a(3,3)

    determinant = a(1,1)*(a(3,3)*a(2,2)-a(3,2)*a(2,3)) &
      - a(2,1)*(a(3,3)*a(1,2)-a(3,2)*a(1,3)) &
      + a(3,1)*(a(2,3)*a(1,2)-a(2,2)*a(1,3))
  end function determinant

  !----------------------------------------------------
  pure function invert_cell(a) result(b)
    real(kind=dp), intent(in)  :: a(3,3)
    real(kind=dp)  :: b(3,3)

    real(kind=dp) :: idet

    idet=1.0_dp/determinant(a)
    b(1,1)= idet*(a(3,3)*a(2,2)-a(3,2)*a(2,3))
    b(1,2)= idet*(a(3,1)*a(2,3)-a(3,3)*a(2,1))
    b(1,3)= idet*(a(3,2)*a(2,1)-a(3,1)*a(2,2))

    b(2,1)= idet*(a(3,2)*a(1,3)-a(3,3)*a(1,2))
    b(2,2)= idet*(a(3,3)*a(1,1)-a(3,1)*a(1,3))
    b(2,3)= idet*(a(3,1)*a(1,2)-a(3,2)*a(1,1))

    b(3,1)= idet*(a(2,3)*a(1,2)-a(2,2)*a(1,3))
    b(3,2)= idet*(a(2,1)*a(1,3)-a(2,3)*a(1,1))
    b(3,3)= idet*(a(2,2)*a(1,1)-a(2,1)*a(1,2))
  end function invert_cell

  ! init random numbers, optionally with a seed
  ! to generate them just use call randomn_number(d)
  ! with d scalar or array
  subroutine init_random()

    integer(kind=ip), allocatable :: seed(:)
    integer(kind=ip) :: p,is,i

    is=fdf_integer('seed', 13)
    call random_seed(size=p)
    allocate(seed(p))
    seed=17*[(i-is,i=1,p)]
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine init_random

end module numeric_esl
