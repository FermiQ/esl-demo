module esl_grid_m
  use prec, only : dp,lp
  use iso_c_binding

  use esl_numeric_m, only: grylmr
  use module_fft_sg
  use pspiof_m
  
  implicit none
  include 'fftw3.f03'

  private

  public :: grid_t,      &
            integrate,   &
            rs_cube2grid,&
            rs_grid2cube,&
            overlap

  !Data structure for the real space grid
  type grid_t
    real(dp) :: hgrid(3) !< Real space spacing
    integer :: ndims(3)  !< Number of points in each directions
    integer :: np !< Total number of points in the real space grid
    real(dp), allocatable :: r(:,:) !<Grid point coordinates 
    real(dp) :: volelem !<Volume element
    real(dp) :: volume

    type(C_PTR) :: fftplan !< Forward FFT plan
    type(C_PTR) :: ifftplan !< Backward FFT (IFFT) plan
  contains
    
    private
    procedure, public :: init
    procedure, public :: radial_function
    procedure, public :: radial_function_gradient
    procedure, public :: radial_function_ylm_gradient
    procedure, public :: radial_function_ylm
    procedure, public :: overlap
    procedure, public :: matrix_elem
    procedure, public :: summary
    
    procedure, private :: dintegrate, zintegrate
    generic, public :: integrate => dintegrate, zintegrate
    
    final :: cleanup
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
    complex(dp), allocatable :: rin(:,:,:), cout(:,:,:)

    this%ndims = ndims

    do idim = 1, 3
      do ! Ensure that Poisson Solver double grid will work.
        call fourier_dim(2 * this%ndims(idim), twice)
        if (2 * this%ndims(idim) == twice) exit
        call fourier_dim(this%ndims(idim) + 1, this%ndims(idim))
      end do
      this%hgrid(idim) = cell(idim, idim) / real(this%ndims(idim)+1, dp)
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
    this%volume = cell(1,1)*(cell(2,2)*cell(3,3)-cell(2,3)*cell(3,2)) - &
          cell(1,2)*(cell(2,1)*cell(3,3)-cell(2,3)*cell(3,1)) + &
          cell(1,3)*(cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1))

    ! Initialization for FFT and IFFT
    allocate(rin(this%ndims(1), this%ndims(2), this%ndims(3)))
    allocate(cout(this%ndims(1), this%ndims(2), this%ndims(3)))
    this%fftplan = fftw_plan_dft_3d(this%ndims(1), this%ndims(2), this%ndims(3), &
      rin, cout, FFTW_FORWARD, FFTW_ESTIMATE)
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
    call yaml_map('Volume element', this%volelem)
    call yaml_mapping_close()

  end subroutine summary

  ! Evaluate a radial function on the real-space grid
  !----------------------------------------------------
  subroutine radial_function(this, rfunc, r_center, func)

    class(grid_t) :: this
    type(pspiof_meshfunc_t), intent(in)  :: rfunc
    real(dp),                intent(in)  :: r_center(3)
    real(dp),                intent(out) :: func(:)

    integer :: ip
    real(dp) :: xyz(3), r

    do ip = 1, this%np
        
      xyz(:) = this%r(:,ip) - r_center(:)
      r = sqrt(sum(xyz**2))
      func(ip) = pspiof_meshfunc_eval(rfunc, r)
        
    end do

  end subroutine radial_function

  ! Evaluate the gradient of a radial function on the real-space grid
  !----------------------------------------------------
  subroutine radial_function_gradient(this, rfunc, r_center, gfunc)
    class(grid_t) :: this
    type(pspiof_meshfunc_t), intent(in)  :: rfunc
    real(dp),                intent(in)  :: r_center(3)
    real(dp),                intent(out) :: gfunc(:,:)

    integer :: ip
    real(dp) :: xyz(3), r
      
    do ip = 1, this%np

      xyz(:) = this%r(:,ip) - r_center(:)
      r = sqrt(sum(xyz**2))
      gfunc(1:3, ip) = pspiof_meshfunc_eval_deriv(rfunc, r)
        
    end do

  end subroutine radial_function_gradient

  ! Evaluate a radial function times a spherical harmonic on the real-space grid
  !----------------------------------------------------
  subroutine radial_function_ylm(this, rfunc, ll, mm, r_center, func)
    class(grid_t) :: this
    type(pspiof_meshfunc_t), intent(in)  :: rfunc
    integer,                 intent(in)  :: ll
    integer,                 intent(in)  :: mm
    real(dp),                intent(in)  :: r_center(3)
    real(dp),                intent(out) :: func(:)

    integer :: ip
    real(dp) :: xyz(3), r
      
    do ip = 1, this%np
      
      xyz(:) = this%r(:,ip) - r_center(:)
      call grylmr(xyz(1), xyz(2), xyz(3), ll, mm, func(ip))

      r = sqrt(sum(xyz**2))
      func(ip) = func(ip)*pspiof_meshfunc_eval(rfunc, r)
        
    end do
     
  end subroutine radial_function_ylm

  ! Evaluate the gradient of a radial function times a spherical harmonic on the real-space grid
  !----------------------------------------------------
  subroutine radial_function_ylm_gradient(this, rfunc, ll, mm, r_center, gfunc)
    class(grid_t) :: this
    type(pspiof_meshfunc_t), intent(in)  :: rfunc
    integer,                 intent(in)  :: ll
    integer,                 intent(in)  :: mm
    real(dp),                intent(in)  :: r_center(3)
    real(dp),                intent(out) :: gfunc(:,:)

    integer :: ip
    real(dp) :: xyz(3), r, f

    do ip = 1, this%np

      xyz(:) = this%r(:,ip) - r_center(:)
      call grylmr(xyz(1), xyz(2), xyz(3), ll, mm, f, gfunc(1:3,ip)) 

      r = sqrt(sum(xyz**2))
      gfunc(1:3, ip) = f*pspiof_meshfunc_eval_deriv(rfunc, r) + &
          gfunc(1:3, ip)*pspiof_meshfunc_eval(rfunc, r)
        
    end do
     
  end subroutine radial_function_ylm_gradient

  !Integrate a function over the real-space grid
  !----------------------------------------------------
  function dintegrate(grid, ff) result(int)
    class(grid_t), intent(in) :: grid
    real(dp), intent(in) :: ff(:)
    real(dp) :: int

    integer :: ip

    int = 0._dp
    do ip = 1 , grid%np
      int = int + ff(ip)
    end do
    int = int*grid%volelem
    
  end function dintegrate

  !Integrate a function over a complex-space grid
  !----------------------------------------------------
  function zintegrate(grid, ff) result(int)
    class(grid_t), intent(in) :: grid
    complex(dp), intent(in) :: ff(:)
    complex(dp) :: int

    integer :: ip

    int = cmplx(0._dp, 0._dp, dp)
    do ip = 1 , grid%np
       int = int + ff(ip)
    end do
    int = int*grid%volelem
    
  end function zintegrate


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
    overlap = 0._dp
    if ( dist < r1 + r2 ) then

      do ip = 1 , grid%np
        overlap = overlap + ao1(ip)*ao2(ip)
      end do
      overlap = overlap*grid%volelem

    end if

  end function overlap

  !Matrix element <AO-1 | pot | AO-2 >
  !----------------------------------------------------
  function matrix_elem(grid, xyz1, ao1, r1, pot, xyz2, ao2, r2) result(M)
    class(grid_t), intent(in) :: grid
    real(dp), intent(in) :: xyz1(3)
    real(dp), intent(in) :: ao1(:)
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: pot(:)
    real(dp), intent(in) :: xyz2(3)
    real(dp), intent(in) :: ao2(:)
    real(dp), intent(in) :: r2
    real(dp) :: M

    integer :: ip
    real(dp) :: dist

    dist = sqrt(sum( (xyz1 - xyz2) ** 2 ))
    M = 0._dp
    if ( dist < r1 + r2 ) then

       do ip = 1 , grid%np
         M = M + ao1(ip) * pot(ip) * ao2(ip)
       end do
       M = M*grid%volelem
       
    end if

  end function matrix_elem


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
