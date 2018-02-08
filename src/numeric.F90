module numeric_esl
  use prec, only : dp,ip
  use fdf, only : fdf_integer

  implicit none
  private

  public ::              &
            matr3inv,    &
            init_random, &
            grylmr,      &
            distance

contains

  !----------------------------------------------------
  pure real(kind=dp) function determinant(a)
    real(kind=dp), intent(in)  :: a(3,3)

    determinant = a(1,1)*(a(3,3)*a(2,2)-a(3,2)*a(2,3)) &
      - a(2,1)*(a(3,3)*a(1,2)-a(3,2)*a(1,3)) &
      + a(3,1)*(a(2,3)*a(1,2)-a(2,2)*a(1,3))
  end function determinant

  !----------------------------------------------------
  pure function matr3inv(a) result(b)
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
  end function matr3inv

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


  ! ---------------------------------------------------------
  !> This is a Numerical Recipes-based subroutine
  !! computes real spherical harmonics ylm at position (x, y, z)
  !! and the gradients of ylm:
  !!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
  !!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
  !! with (theta,phi) the polar angles of r, c a positive normalization
  subroutine grylmr(x, y, z, li, mi, ylm, grylm)
    use numerics, only : pi
    real(kind=dp),           intent(in)  :: x, y, z
    integer,                 intent(in)  :: li, mi
    real(kind=dp),           intent(out) :: ylm
    real(kind=dp), optional, intent(out) :: grylm(3)

    integer, parameter :: lmaxd = 20
    real(kind=dp),   parameter :: tiny = 1.e-30
    integer :: i, ilm0, l, m, mabs
    integer, save :: lmax = -1

    real(kind=dp) :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac, &
      fourpi, plgndr, phase, pll, pmm, pmmp1, sinm, &
      sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize
    real(kind=dp), save :: c(0:(lmaxd+1)*(lmaxd+1))

    ! evaluate normalization constants once and for all
    if (li > lmax) then
      fourpi = 4.d0*PI
      do l = 0, li
        ilm0 = l*l + l
        do m = 0, l
          fac = (2*l+1)/fourpi
          do i = l - m + 1, l + m
            fac = fac/i
          end do
          c(ilm0 + m) = sqrt(fac)
          ! next line because ylm are real combinations of m and -m
          if(m /= 0) c(ilm0 + m) = c(ilm0 + m)*sqrt(2.d0)
          c(ilm0 - m) = c(ilm0 + m)
        end do
      end do
      lmax = li
    end if

    ! if l=0, no calculations are required
    if (li == 0) then
      ylm = c(0)
      if(present(grylm)) grylm(:) = 0.d0
      return
    end if

    ! if r=0, direction is undefined => make ylm=0 except for l=0
    r2 = x**2 + y**2 + z**2
    if(r2 < tiny) then
      ylm = 0.d0
      if(present(grylm)) grylm(:) = 0.d0
      return
    end if
    rsize = sqrt(r2)

    Rx = x/rsize
    Ry = y/rsize
    Rz = z/rsize

    ! explicit formulas for l=1 and l=2
    if(li == 1) then
      select case(mi)
      case(-1)
        ylm = (-c(1))*Ry
        if(present(grylm)) then
          grylm(1) = c(1)*Rx*Ry/rsize
          grylm(2) = (-c(1))*(1.d0 - Ry*Ry)/rsize
          grylm(3) = c(1)*Rz*Ry/rsize
        end if
      case(0)
        ylm = c(2)*Rz
        if(present(grylm)) then
          grylm(1) = (-c(2))*Rx*Rz/rsize
          grylm(2) = (-c(2))*Ry*Rz/rsize
          grylm(3) = c(2)*(1.d0 - Rz*Rz)/rsize
        end if
      case(1)
        ylm = (-c(3))*Rx
        if(present(grylm)) then
          grylm(1) = (-c(3))*(1.d0 - Rx*Rx)/rsize
          grylm(2) = c(3)*Ry*Rx/rsize
          grylm(3) = c(3)*Rz*Rx/rsize
        end if
      end select
      return
    end if

    if(li == 2) then
      select case(mi)
      case(-2)
        ylm = c(4)*6.d0*Rx*Ry
        if(present(grylm)) then
          grylm(1) = (-c(4))*6.d0*(2.d0*Rx*Rx*Ry - Ry)/rsize
          grylm(2) = (-c(4))*6.d0*(2.d0*Ry*Rx*Ry - Rx)/rsize
          grylm(3) = (-c(4))*6.d0*(2.d0*Rz*Rx*Ry)/rsize
        end if
      case(-1)
        ylm = (-c(5))*3.d0*Ry*Rz
        if(present(grylm)) then
          grylm(1) = c(5)*3.d0*(2.d0*Rx*Ry*Rz)/rsize
          grylm(2) = c(5)*3.d0*(2.d0*Ry*Ry*Rz - Rz)/rsize
          grylm(3) = c(5)*3.d0*(2.d0*Rz*Ry*Rz - Ry)/rsize
        end if
      case(0)
        ylm = c(6)*0.5d0*(3.d0*Rz*Rz - 1.d0)
        if(present(grylm)) then
          grylm(1) = (-c(6))*3.d0*(Rx*Rz*Rz)/rsize
          grylm(2) = (-c(6))*3.d0*(Ry*Rz*Rz)/rsize
          grylm(3) = (-c(6))*3.d0*(Rz*Rz - 1.d0)*Rz/rsize
        end if
      case(1)
        ylm = (-c(7))*3.d0*Rx*Rz
        if(present(grylm)) then
          grylm(1) = c(7)*3.d0*(2.d0*Rx*Rx*Rz - Rz)/rsize
          grylm(2) = c(7)*3.d0*(2.d0*Ry*Rx*Rz)/rsize
          grylm(3) = c(7)*3.d0*(2.d0*Rz*Rx*Rz - Rx)/rsize
        end if
      case(2)
        ylm = c(8)*3.d0*(Rx*Rx - Ry*Ry)
        if(present(grylm)) then
          grylm(1) = (-c(8))*6.d0*(Rx*Rx - Ry*Ry - 1.d0)*Rx/rsize
          grylm(2) = (-c(8))*6.d0*(Rx*Rx - Ry*Ry + 1.d0)*Ry/rsize
          grylm(3) = (-c(8))*6.d0*(Rx*Rx - Ry*Ry)*Rz/rsize
        end if
      end select
      return
    end if

    ! general algorithm based on routine plgndr of numerical recipes
    mabs = abs(mi)
    xysize = sqrt(max(Rx*Rx + Ry*Ry, tiny))
    cosphi = Rx/xysize
    sinphi = Ry/xysize
    cosm = 1.d0
    sinm = 0.d0
    do m = 1, mabs
      cosmm1 = cosm
      sinmm1 = sinm
      cosm = cosmm1*cosphi - sinmm1*sinphi
      sinm = cosmm1*sinphi + sinmm1*cosphi
    end do

    if(mi < 0) then
      phase = sinm
      dphase = mabs*cosm
    else
      phase = cosm
      dphase = (-mabs)*sinm
    end if

    pmm = 1.d0
    fac = 1.d0

    if(mabs > 0.d0) then
      do i = 1, mabs
        pmm = (-pmm)*fac*xysize
        fac = fac + 2.d0
      end do
    end if

    if(li == mabs) then
      plgndr = pmm
      dplg = (-li)*Rz*pmm/(xysize**2)
    else
      pmmp1 = Rz*(2*mabs + 1)*pmm
      if(li == mabs + 1) then
        plgndr = pmmp1
        dplg = -((li*Rz*pmmp1 - (mabs + li)*pmm)/(xysize**2))
      else
        do l = mabs + 2, li
          pll = (Rz*(2*l - 1)*pmmp1 - (l + mabs - 1)*pmm)/(l - mabs)
          pmm = pmmp1
          pmmp1 = pll
        end do
        plgndr = pll
        dplg = -((li*Rz*pll - (l + mabs - 1)*pmm)/(xysize**2))
      end if
    end if

    ilm0 = li*li + li
    cmi = c(ilm0 + mi)
    ylm = cmi*plgndr*phase

    if(present(grylm)) then
      grylm(1) = (-cmi)*dplg*Rx*Rz*phase/rsize     &
        -cmi*plgndr*dphase*Ry/(rsize*xysize**2)
      grylm(2) = (-cmi)*dplg*Ry*Rz*phase/rsize     &
        +cmi*plgndr*dphase*Rx/(rsize*xysize**2)
      grylm(3)= cmi*dplg*(1.d0 - Rz*Rz)*phase/rsize
    end if

    return
  end subroutine grylmr

  real(kind=dp) function distance(r1, r2)
    real(kind=dp), intent(in) :: r1(3)
    real(kind=dp), intent(in) :: r2(3)

    distance = sqrt(sum((r1(1:3)-r2(1:3))**2))

  end function distance

end module numeric_esl
