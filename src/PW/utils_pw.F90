module esl_utils_pw_m

  use prec
  use yaml_output
  use esl_geometry_m
  use esl_grid_m
  use iso_c_binding

  implicit none
  include 'fftw3.f03'
  
  private

  public ::                          &
            get_number_of_pw,        &
            pw2grid,                 &
            construct_mod_map_tables,&
            grid2pw 
            

contains


  integer function get_number_of_pw(ndims, ecut, gmet, kpt) result(npw)
    use esl_constants_m, only: pi

    integer,          intent(in) :: ndims(3)
    real(dp),         intent(in) :: ecut
    real(dp),         intent(in) :: gmet(3, 3)
    real(dp),         intent(in) :: kpt(3)

    integer :: i1, i2, i3
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

  subroutine construct_mod_map_tables(ndims, ecut, gmet, kpt, gmod2, gmap)
    use esl_constants_m, only: pi

    integer,          intent(in) :: ndims(3)
    real(dp),         intent(in) :: ecut
    real(dp),         intent(in) :: gmet(3, 3)
    real(dp),         intent(in) :: kpt(3)
    real(dp),        intent(out) :: gmod2(:)
    integer,         intent(out) :: gmap(:,:)

    integer :: i1, i2, i3, npw
    real(dp) :: threshold, tmp_gmod

    npw = 0
    threshold = 0.5_dp * ecut / pi**2
    do i1 = -ndims(1) / 2, ndims(1)/2
       do i2 = -ndims(2) / 2, ndims(2)/2
          do i3 = -ndims(3) / 2, ndims(3)/2
             tmp_gmod = dsq(gmet, kpt, i1, i2, i3)
             if (tmp_gmod <= threshold) then
               npw = npw + 1
               gmod2(npw) = tmp_gmod
               gmap(1,npw) = i1 + ndims(1)/2 + 1
               gmap(2,npw) = i2 + ndims(2)/2 + 1
               gmap(3,npw) = i3 + ndims(3)/2 + 1
             endif
          end do
       end do
    end do

  end subroutine construct_mod_map_tables
  
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
    

  subroutine pw2grid(grid, gmap, ndims, npw, coef_pw, coef_rs)
    type(grid_t),          intent(in) :: grid
    integer,               intent(in) :: gmap(:,:) 
    integer,               intent(in) :: ndims(3)
    integer,               intent(in) :: npw
    complex(kind=dp),      intent(in) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),     intent(out) :: coef_rs(:) !(np)

    complex(kind=dp), allocatable :: fourier_cube(:,:,:)
    complex(kind=dp), allocatable :: rs_cube(:,:,:)
    
    allocate(fourier_cube(1:ndims(1),1:ndims(2),1:ndims(3)))
    allocate(rs_cube(1:ndims(1),1:ndims(2),1:ndims(3)))

    call fourier_sphere2cube(gmap, ndims, npw, coef_pw, fourier_cube)    

    ! FFT-1
    call fftw_execute_dft(grid%iFFTplan, fourier_cube(:,:,:), rs_cube(:,:,:))

    call rs_cube2grid(grid, rs_cube, coef_rs)

    deallocate(fourier_cube, rs_cube)

  end subroutine pw2grid

  subroutine grid2pw(grid, gmap, ndims, npw, coef_rs, coef_pw, accumulate)
    type(grid_t),          intent(in) :: grid
    integer,               intent(in) :: gmap(:,:)
    integer,               intent(in) :: ndims(3)
    integer,               intent(in) :: npw
    complex(kind=dp),   intent(inout) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),      intent(in) :: coef_rs(:) !(np)
    logical, optional,     intent(in) :: accumulate

    complex(kind=dp), allocatable :: fourier_cube(:,:,:)
    complex(kind=dp), allocatable :: rs_cube(:,:,:)

    allocate(fourier_cube(1:ndims(1),1:ndims(2),1:ndims(3)))
    allocate(rs_cube(1:ndims(1),1:ndims(2),1:ndims(3)))

    call rs_grid2cube(grid, coef_rs, rs_cube)

    ! FFT
    call fftw_execute_dft(grid%FFTplan, rs_cube(:,:,:), fourier_cube(:,:,:))

    call fourier_cube2sphere(gmap, npw, fourier_cube, coef_pw, accumulate)

    deallocate(fourier_cube, rs_cube)

  end subroutine grid2pw


  subroutine fourier_sphere2cube(gmap, ndims, npw, coef_pw, fourier_cube)
    integer,               intent(in) :: gmap(:,:)
    integer,               intent(in) :: ndims(3)
    integer,               intent(in) :: npw
    complex(kind=dp),      intent(in) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),     intent(out) :: fourier_cube(:,:,:)

    integer :: ipw

    fourier_cube(1:ndims(1), 1:ndims(2), 1:ndims(3)) = cmplx(0.d0,0.d0,kind=dp)

    do ipw = 1, npw
      fourier_cube(gmap(1,ipw), gmap(2,ipw), gmap(3,ipw)) = coef_pw(ipw)
    end do    

  end subroutine fourier_sphere2cube

  subroutine fourier_cube2sphere(gmap, npw, fourier_cube, coef_pw, accumulate)
    integer,               intent(in) :: gmap(:,:)
    integer,               intent(in) :: npw
    complex(kind=dp),   intent(inout) :: coef_pw(:) !(pw%npw)
    complex(kind=dp),      intent(in) :: fourier_cube(:,:,:)
    logical, optional,     intent(in) :: accumulate

    integer :: ipw
    logical :: acc

    acc = .false.
    if(present(accumulate)) acc = accumulate 

    if(.not.acc) then
      coef_pw(1:npw) = cmplx(0.d0,0.d0,kind=dp)
 
      do ipw = 1, npw
        coef_pw(ipw) = fourier_cube(gmap(1,ipw), gmap(2,ipw), gmap(3,ipw))
      end do
    else
      do ipw = 1, npw
        coef_pw(ipw) = coef_pw(ipw) + fourier_cube(gmap(1,ipw), gmap(2,ipw), gmap(3,ipw))
      end do
    end if

  end subroutine fourier_cube2sphere


end module esl_utils_pw_m
