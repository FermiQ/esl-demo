module esl_grid_m
  use prec, only : dp,ip
  use module_fft_sg
  implicit none

  private

  public :: grid_t,   &
            integrate

  !Data structure for the real space grid
  type grid_t
    real(dp) :: hgrid(3) !< Real space spacing
    integer  :: ndims(3)  !< Number of points in each directions
    integer  :: np !< Total number of points in the real space grid
    real(dp), allocatable :: r(:,:) !<Grid point coordinates 
    real(dp) :: volelem !<Volume element
  contains
    private
    procedure, public :: init
    procedure, public :: get_atomic_orbital
    procedure, public :: summary
    final  :: cleanup
  end type grid_t

  interface integrate
    module procedure dintegrate, zintegrate
  end interface integrate

contains

  !Initialize the grid
  !----------------------------------------------------
  subroutine init(this, ndims, cell)
    class(grid_t) :: this
    integer,  intent(inout) :: ndims(3)
    real(dp), intent(in) :: cell(3,3)

    integer :: idim, ix, iy, iz, ip
    integer :: n, twice

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

    ndims = this%ndims
  end subroutine init

  !Release the grid
  !----------------------------------------------------
  subroutine cleanup(this)
    type(grid_t) :: this

    if(allocated(this%r)) deallocate(this%r)

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
  subroutine get_atomic_orbital(this, ll, mm, r_at, ao, grad_ao)
    use esl_numeric_m, only: grylmr
    class(grid_t) :: this
    integer,        intent(in) :: ll
    integer,        intent(in) :: mm
    real(dp), intent(out) :: r_at(3)
    real(dp), intent(out) :: ao(:)
    real(dp), intent(out) :: grad_ao(:,:)

    integer :: ip
    real(dp) :: x, y, z, r

    do ip = 1, this%np
       x = this%r(1,ip) - r_at(1)
       y = this%r(2,ip) - r_at(2)
       z = this%r(3,ip) - r_at(3)
       call grylmr(x, y, z, ll, mm, ao(ip), grad_ao(1:3,ip)) 

       r = sqrt(x**2+y**2+z**2)
       !Here we need to multiply by the radial part
    end do

  end subroutine get_atomic_orbital

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

    int_ff = cmplx(0.d0,0.d0)
    do ip = 1,grid%np
       int_ff = int_ff + ff(ip)
    end do
    int_ff = int_ff*grid%volelem

  end subroutine zintegrate


  !Overlap
  !----------------------------------------------------
  real(dp) function overlap(grid, xyz1, ao1, r1, xyz2, ao2, r2)
    type(grid_t),    intent(in) :: grid
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

end module esl_grid_m
