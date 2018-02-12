module esl_utils_pw_m

  use prec
  use yaml_output
  use esl_geometry_m
  use esl_grid_m
  
  private

  public ::                          &
            get_number_of_pw,        &
            pw2grid,                 &
            grid2pw 
            

contains


  integer function get_number_of_pw(ndims, ecut, gmet, kpt) result(npw)
    use esl_constants_m, only: pi

    integer,          intent(in) :: ndims(3)
    real(dp),         intent(in) :: ecut
    real(dp),         intent(in) :: gmet(3, 3)
    real(dp),         intent(in) :: kpt(3)

    integer :: i1, i2, i3, i
    real(dp) :: threshold

    npw = 0
    threshold = 0.5_dp * ecut / pi**2
    do i1 = -ndims(1) / 2, ndims(1)/2
       do i2 = -ndims(2) / 2, ndims(2)/2
          do i3 = -ndims(3) / 2, ndims(3)/2
             if (dsq(gmet, kpt, i1, i2, i3) <= threshold) npw = npw + 1
          end do
       end do
    end do

  end function get_number_of_pw  
  
  function dsq(gmet, kpt, i1, i2, i3)
      real(kind=dp), intent(in) :: gmet(3,3)
      real(kind=dp), intent(in) :: kpt(3)
      integer,       intent(in) :: i1, i2, i3
      real(dp) :: dsq

      dsq = gmet(1, 1)*(kpt(1) + dble(i1))**2 &
           & + gmet(2, 2)*(kpt(2) + dble(i2))**2 &
           & + gmet(3, 3)*(kpt(3) + dble(i3))**2 &
           & + 2._dp*(gmet(1, 2)*(kpt(1) + dble(i1))*(kpt(2) + dble(i2)) &
           & + gmet(2, 3)*(kpt(2) + dble(i2))*(kpt(3) + dble(i3)) &
           & + gmet(3, 1)*(kpt(3) + dble(i3))*(kpt(1) + dble(i1)))

  end function dsq
    

  subroutine pw2grid(grid, gmet, kpt, ndims, ecut, coef_pw, coef_rs)
    type(grid_t),          intent(in) :: grid
    real(kind=dp),         intent(in) :: gmet(3,3)
    real(kind=dp),         intent(in) :: kpt(3) 
    integer,               intent(in) :: ndims(3)
    real(kind=dp),         intent(in) :: ecut
    complex(kind=dp),      intent(in) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),     intent(out) :: coef_rs(:) !(np)

    complex(kind=dp), allocatable :: fourier_cube(:,:,:)
    complex(kind=dp), allocatable :: rs_cube(:,:,:)
    
    allocate(fourier_cube(1:ndims(1),1:ndims(2),1:ndims(3)))
    allocate(rs_cube(1:ndims(1),1:ndims(2),1:ndims(3)))

    call fourier_sphere2cube(gmet, kpt, ndims, ecut, coef_pw, fourier_cube)    

    ! FFT-1
    ! TODO include 'fftw3.f90' should be put properly
    call dfftw_execute_dft(grid%iFFTplan, fourier_cube, fourier_cube)

    ! FFT
    call dfftw_execute_dft(grid%FFTplan, fourier_cube, fourier_cube)

    call rs_cube2grid(grid, rs_cube, coef_rs)

    deallocate(fourier_cube, rs_cube)

  end subroutine pw2grid

  subroutine grid2pw(grid, gmet, kpt, ndims, ecut, coef_rs, coef_pw)
    type(grid_t),          intent(in) :: grid
    real(kind=dp),         intent(in) :: gmet(3,3)
    real(kind=dp),         intent(in) :: kpt(3)
    integer,               intent(in) :: ndims(3)
    real(kind=dp),         intent(in) :: ecut
    complex(kind=dp),     intent(out) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),      intent(in) :: coef_rs(:) !(np)

    complex(kind=dp), allocatable :: fourier_cube(:,:,:)
    complex(kind=dp), allocatable :: rs_cube(:,:,:)

    allocate(fourier_cube(1:ndims(1),1:ndims(2),1:ndims(3)))
    allocate(rs_cube(1:ndims(1),1:ndims(2),1:ndims(3)))

    call rs_grid2cube(grid, coef_rs, rs_cube)

    !Here FFT

    call fourier_cube2sphere(gmet, kpt, ndims, ecut, fourier_cube, coef_pw)

    deallocate(fourier_cube, rs_cube)

  end subroutine grid2pw


  subroutine fourier_sphere2cube(gmet, kpt, ndims, ecut, coef_pw, fourier_cube)
    use esl_constants_m, only: pi
    real(kind=dp),         intent(in) :: gmet(3,3)
    real(kind=dp),         intent(in) :: kpt(3)
    integer,               intent(in) :: ndims(3)
    real(kind=dp),         intent(in) :: ecut
    complex(kind=dp),      intent(in) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),     intent(out) :: fourier_cube(:,:,:)

    integer :: i1, i2, i3, npw
    real(dp) :: threshold

    fourier_cube(1:ndims(1), 1:ndims(2), 1:ndims(3)) = 0.d0

    npw = 0
    threshold = 0.5_dp * ecut / pi**2
    do i1 = 1, ndims(1)
       do i2 = 1, ndims(2)
          do i3 = 1, ndims(3)
            if (dsq(gmet, kpt, i1-ndims(1)/2-1, i2-ndims(2)/2-1, i3-ndims(3)/2-1) <= threshold) then
               npw = npw + 1
               fourier_cube(i1, i2, i3) = coef_pw(npw)
            end if
          end do
       end do
    end do    

  end subroutine fourier_sphere2cube

  subroutine fourier_cube2sphere(gmet, kpt, ndims, ecut, fourier_cube, coef_pw)
    use esl_constants_m, only: pi  
    real(kind=dp),         intent(in) :: gmet(3,3)
    real(kind=dp),         intent(in) :: kpt(3)
    integer,               intent(in) :: ndims(3)
    real(kind=dp),         intent(in) :: ecut
    complex(kind=dp),     intent(out) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),      intent(in) :: fourier_cube(:,:,:)

    integer :: i1, i2, i3, npw
    real(dp) :: threshold

    coef_pw(:) = 0.d0

    npw = 0
    threshold = 0.5_dp * ecut / pi**2
    do i1 = 1, ndims(1)
       do i2 = 1, ndims(2)
          do i3 = 1, ndims(3)
            if (dsq(gmet, kpt, i1-ndims(1)/2-1, i2-ndims(2)/2-1, i3-ndims(3)/2-1) <= threshold) then
               npw = npw + 1
               coef_pw(npw) = fourier_cube(i1, i2, i3)
            end if
          end do
       end do
    end do

  end subroutine fourier_cube2sphere


end module esl_utils_pw_m
