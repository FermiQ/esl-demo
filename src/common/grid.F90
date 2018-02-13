module esl_grid_m
  use prec, only : dp,lp
  use iso_c_binding
  use module_fft_sg
  use pspiof_m
  implicit none
  include 'fftw3.f03'

  private

  public :: grid_t,      &
            integrate,   &
            rs_cube2grid,&
            rs_grid2cube,&
            overlap,     &
            matrixelem

  !Data structure for the real space grid
  type grid_t
    real(dp) :: hgrid(3) !< Real space spacing
    integer  :: ndims(3)  !< Number of points in each directions
    integer  :: np !< Total number of points in the real space grid
    real(dp), allocatable :: r(:,:) !<Grid point coordinates 
    real(dp) :: volelem !<Volume element

    type(C_PTR) fftplan !< Forward FFT plan
    type(C_PTR) ifftplan !< Backward FFT (IFFT) plan
  contains
    private
    procedure, public :: init
    procedure, public :: radial_function
    procedure, public :: overlap
    procedure, public :: summary
    final  :: cleanup
  end type grid_t

  interface integrate
    module procedure dintegrate, zintegrate
  end interface integrate

  interface rs_cube2grid
    module procedure drs_cube2grid, zrs_cube2grid
  end interface rs_cube2grid

  interface rs_grid2cube
    module procedure drs_grid2cube, zrs_grid2cube
  end interface rs_grid2cube


contains

  !Initialize the grid
  !----------------------------------------------------
  subroutine init(this, ndims, cell)
    class(grid_t) :: this
    integer,  intent(inout) :: ndims(3)
    real(dp), intent(in) :: cell(3,3)

    integer :: idim, ix, iy, iz, ip
    integer :: twice
    complex(dp),         allocatable :: rin(:,:,:), cout(:,:,:)

    this%ndims = ndims

    do idim = 1, 3
      do ! Ensure that Poisson Solver double grid will work.
        call fourier_dim(2 * this%ndims(idim), twice)
        if (2 * this%ndims(idim) == twice) exit
        call fourier_dim(this%ndims(idim) + 1, this%ndims(idim))
      end do
      this%hgrid(idim) = cell(idim, idim) / real(this%ndims(idim) + 1, dp)
    end do
    this%np = this%ndims(1)*this%ndims(2)*this%ndims(3)

    !Generation of the grid points
    allocate(this%r(3, this%np))
    ip = 0
    do ix = 1, this%ndims(1)
      do iy = 1, this%ndims(2)
        do iz = 1, this%ndims(3)
          ip = ip + 1
          this%r(1, ip) = ix*this%hgrid(1) - 0.5d0*cell(1,1)
          this%r(2, ip) = iy*this%hgrid(2) - 0.5d0*cell(2,2)
          this%r(3, ip) = iz*this%hgrid(3) - 0.5d0*cell(3,3)
        end do
      end do
    end do
    !We have a cubic cell
    this%volelem = this%hgrid(1)*this%hgrid(2)*this%hgrid(3)

    ! Initialization for FFT and IFFT
    allocate(rin(this%ndims(1), this%ndims(2), this%ndims(3)))
    allocate(cout(this%ndims(1)/2+1, this%ndims(2), this%ndims(3)))
    this%fftplan = fftw_plan_dft_3d(this%ndims(1), this%ndims(2), this%ndims(3), &
      rin, cout, FFTW_FORWARD, FFTW_ESTIMATE)
    deallocate(rin,cout)

    allocate(rin(this%ndims(1)/2+1, this%ndims(2), this%ndims(3)))
    allocate(cout(this%ndims(1), this%ndims(2), this%ndims(3)))
    this%ifftplan = fftw_plan_dft_3d(this%ndims(1), this%ndims(2), this%ndims(3), &
      rin, cout, FFTW_BACKWARD, FFTW_ESTIMATE)
    deallocate(rin,cout)

    ndims = this%ndims
  end subroutine init

  !Release the grid
  !----------------------------------------------------
  subroutine cleanup(this)
    type(grid_t) :: this

    if(allocated(this%r)) deallocate(this%r)

    ! Deconstructor for fft plan
    call dfftw_destroy_plan(this%fftplan)
    call dfftw_destroy_plan(this%ifftplan)

  end subroutine cleanup

  !summary
  !----------------------------------------------------
  subroutine summary(this)
    use yaml_output
    class(grid_t) :: this

    call yaml_mapping_open("Grid")
    call yaml_map("Spacing", this%hgrid)
    call yaml_map("Ndims", this%ndims)
    call yaml_map("Total number of points", this%np)
    call yaml_mapping_close()

  end subroutine summary

  !Evaluate an atomic orbital on the real-space grid
  !----------------------------------------------------
  subroutine radial_function(this, rfunc, ll, mm, r_center, func, gfunc)
    use esl_numeric_m, only: grylmr
    class(grid_t) :: this
    type(pspiof_meshfunc_t), intent(in)  :: rfunc
    integer,                 intent(in)  :: ll
    integer,                 intent(in)  :: mm
    real(dp),                intent(in)  :: r_center(3)
    real(dp),                intent(out) :: func(:)
    real(dp), optional,      intent(out) :: gfunc(:,:)

    integer :: ip
    real(dp) :: x, y, z, r, fr

    do ip = 1, this%np
      x = this%r(1,ip) - r_center(1)
      y = this%r(2,ip) - r_center(2)
      z = this%r(3,ip) - r_center(3)
      if (present(gfunc)) then
        call grylmr(x, y, z, ll, mm, func(ip), gfunc(1:3,ip)) 
      else
        call grylmr(x, y, z, ll, mm, func(ip))
      end if

      r = sqrt(x**2 + y**2 + z**2)
      fr = pspiof_meshfunc_eval(rfunc, r)
      func(ip) = func(ip)*fr
      if (present(gfunc)) then
        gfunc(1:3, ip) = func(ip)*pspiof_meshfunc_eval_deriv(rfunc, r) + gfunc(1:3, ip)*fr
      end if
    end do
     
  end subroutine radial_function

  !Integrate a function over the real-space grid
  !----------------------------------------------------
  subroutine dintegrate(grid, ff, int_ff)
    type(grid_t),    intent(in) :: grid
    real(dp),   intent(in) :: ff(:)
    real(dp),  intent(out) :: int_ff

    integer :: ip

    int_ff = 0.d0
    do ip=1,grid%np
       int_ff = int_ff + ff(ip)
    end do
    int_ff = int_ff*grid%volelem

  end subroutine dintegrate

  !Integrate a function over the real-space grid
  !----------------------------------------------------
  subroutine zintegrate(grid, ff, int_ff)
    type(grid_t),       intent(in) :: grid
    complex(dp),   intent(in) :: ff(:)
    complex(dp),  intent(out) :: int_ff

    integer :: ip

    int_ff = cmplx(0.d0,0.d0, kind=dp)
    do ip = 1,grid%np
       int_ff = int_ff + ff(ip)
    end do
    int_ff = int_ff*grid%volelem

  end subroutine zintegrate


  !Overlap
  !----------------------------------------------------
  real(dp) function overlap(grid, xyz1, ao1, r1, xyz2, ao2, r2)
    class(grid_t),    intent(in) :: grid
    real(dp),   intent(in) :: xyz1(3)
    real(dp),   intent(in) :: ao1(:)
    real(dp),   intent(in) :: r1
    real(dp),   intent(in) :: xyz2(3)
    real(dp),   intent(in) :: ao2(:)
    real(dp),   intent(in) :: r2

    integer :: ip
    real(dp) :: dist
    
    dist = sqrt(sum( (xyz1 - xyz2) ** 2 ))
    if ( dist < r1 + r2 ) then

      do ip = 1 , grid%np
        overlap = overlap + ao1(ip)*ao2(ip)
      end do
      overlap = overlap*grid%volelem

    else
      
      overlap = 0._dp
      
    end if

  end function overlap

  !Matrix element
  !----------------------------------------------------
  real(dp) function matrixelem(grid, xyz1, ao1, r1, xyz2, ao2, r2, pot)
    type(grid_t),    intent(in) :: grid
    real(dp),   intent(in) :: xyz1(3)
    real(dp),   intent(in) :: ao1(:)
    real(dp),   intent(in) :: r1
    real(dp),   intent(in) :: xyz2(3)
    real(dp),   intent(in) :: ao2(:)
    real(dp),   intent(in) :: r2
    real(dp),   intent(in) :: pot(:)

    integer :: ip
    real(dp) :: dist

    dist = sqrt(sum( (xyz1 - xyz2) ** 2 ))
    matrixelem = 0._dp
    if ( dist < r1 + r2 ) then

       do ip = 1 , grid%np
          matrixelem = matrixelem + ao1(ip)*ao2(ip)*pot(ip)
       end do
       matrixelem = matrixelem*grid%volelem

    end if

  end function matrixelem


  subroutine zrs_cube2grid(this, ff_cube, ff_grid)
    type(grid_t),     intent(in)  :: this
    complex(kind=dp), intent(in)  :: ff_cube(:,:,:)
    complex(kind=dp), intent(out) :: ff_grid(:)

    integer :: ip, ix, iy, iz
 
    ip = 0
    do ix = 1, this%ndims(1)
      do iy = 1, this%ndims(2)
        do iz = 1, this%ndims(3)
          ip = ip + 1
          ff_grid(ip) = ff_cube(ix, iy, iz)    
        end do
      end do
    end do

  end subroutine zrs_cube2grid

  subroutine drs_cube2grid(this, ff_cube, ff_grid)
    type(grid_t),  intent(in)  :: this
    real(kind=dp), intent(in)  :: ff_cube(:,:,:)
    real(kind=dp), intent(out) :: ff_grid(:)

    integer :: ip, ix, iy, iz

    ip = 0
    do ix = 1, this%ndims(1)
      do iy = 1, this%ndims(2)
        do iz = 1, this%ndims(3)
          ip = ip + 1
          ff_grid(ip) = ff_cube(ix, iy, iz)
        end do
      end do
    end do

  end subroutine drs_cube2grid

  subroutine zrs_grid2cube(this, ff_grid, ff_cube)
    type(grid_t),  intent(in)  :: this
    complex(kind=dp), intent(out) :: ff_cube(:,:,:)
    complex(kind=dp), intent(in)  :: ff_grid(:)

    integer :: ip, ix, iy, iz

    ip = 0
    do ix = 1, this%ndims(1)
      do iy = 1, this%ndims(2)
        do iz = 1, this%ndims(3)
          ip = ip + 1
          ff_cube(ix, iy, iz) = ff_grid(ip)
        end do
      end do
    end do

  end subroutine zrs_grid2cube


  subroutine drs_grid2cube(this, ff_grid, ff_cube)
    type(grid_t),  intent(in)  :: this
    real(kind=dp), intent(out) :: ff_cube(:,:,:)
    real(kind=dp), intent(in)  :: ff_grid(:)

    integer :: ip, ix, iy, iz

    ip = 0
    do ix = 1, this%ndims(1)
      do iy = 1, this%ndims(2)
        do iz = 1, this%ndims(3)
          ip = ip + 1
          ff_cube(ix, iy, iz) = ff_grid(ip)
        end do
      end do
    end do

  end subroutine drs_grid2cube

end module esl_grid_m
