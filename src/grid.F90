module esl_grid_m
  use prec, only : dp,ip

  use basis_esl

  implicit none

  private

  public :: grid_t
  public :: integrate

  !Data structure for the real space grid
  type grid_t
     real(dp) :: hgrid(3) !< Real space spacing
     integer :: ndims(3)  !< Number of points in each directions
     integer :: np !< Total number of points in the real space grid
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
  subroutine init(this, basis, cell, gcell)
    use module_fft_sg
    class(grid_t) :: this
    type(basis_t), intent(in) :: basis
    real(dp), dimension(3,3), intent(in) :: cell
    real(dp), dimension(3,3), intent(in) :: gcell

    integer :: idim, ix, iy, iz, ip
    integer :: n, twice

    select case (basis%type)
    case (PLANEWAVES)
       call nDimsFromEcut(this%ndims, basis%pw_basis%ecut, gcell, [0._dp, 0._dp, 0._dp])
    case (ATOMICORBS)
       !For the moment the spacing in real space is hardcoded
       !For planewave, this must come from the number of G vectors
       this%hgrid(1:3) = .25
       do idim = 1,3
          n = ceiling(cell(idim, idim) / this%hgrid(idim)) - 1
          call fourier_dim(n, this%ndims(idim))
       end do
    end select

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
    allocate(this%r(3,this%np))
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

  subroutine nDimsFromEcut(ndims, ecut, gcell, kpt)
    use esl_constants_m, only: PI
    use yaml_output
    integer, dimension(3), intent(out) :: ndims
    real(dp), intent(in) :: ecut
    real(dp), dimension(3,3), intent(in) :: gcell
    real(dp), dimension(3), intent(in) :: kpt

    real(dp) :: threshold
    integer :: dir, n, i
    real(dp), dimension(3,3) :: gmet
    real(dp), parameter :: boxcutmin = 2.0_dp

    do i = 1, 3
       gmet(i, :) = gcell(1, i) * gcell(1, :) + &
            gcell(2, i) * gcell(2, :) + &
            gcell(3, i) * gcell(3, :)
    end do

    threshold = 0.5_dp * boxcutmin**2 * ecut / PI**2
    !@todo: Don't take into account symmetries or k points.
    ndims = 16
    do
       if (smallest(ndims, gmet, dir) >= threshold) exit
       call fourier_dim(ndims(dir) + 1, ndims(dir))
    end do

  contains

    function smallest(ndims, gmet, dir)
      integer, dimension(3), intent(in) :: ndims
      real(dp), dimension(3,3), intent(in) :: gmet
      integer, intent(out) :: dir
      real(dp) :: smallest

      integer :: i1, i2, i3, s1, s2, s3, idir
      real(dp) :: prev

      smallest = dsq(ndims(1) / 2, -ndims(2) / 2, -ndims(3) / 2) + 0.01_dp
      do idir = 1, 3
         s1 = 0
         if (idir == 1) s1 = ndims(1)/2
         s2 = 0
         if (idir == 2) s2 = ndims(2)/2
         s3 = 0
         if (idir == 3) s3 = ndims(3)/2
         do i1 = s1, ndims(1)/2
            do i2 = s2, ndims(2)/2
               do i3 = s3, ndims(3)/2
                  prev = smallest
                  smallest = min(smallest, dsq( i1,  i2,  i3))
                  smallest = min(smallest, dsq(-i1,  i2,  i3))
                  smallest = min(smallest, dsq( i1, -i2,  i3))
                  smallest = min(smallest, dsq(-i1, -i2,  i3))
                  smallest = min(smallest, dsq( i1,  i2, -i3))
                  smallest = min(smallest, dsq(-i1,  i2, -i3))
                  smallest = min(smallest, dsq( i1, -i2, -i3))
                  smallest = min(smallest, dsq(-i1, -i2, -i3))
                  if (prev /= smallest) dir = idir
               end do
            end do
         end do
      end do
    end function smallest

    function dsq(i1, i2, i3)
      integer, intent(in) :: i1, i2, i3
      real(dp) :: dsq

      dsq=gmet(1,1)*(kpt(1)+dble(i1))**2&
           & +gmet(2,2)*(kpt(2)+dble(i2))**2&
           & +gmet(3,3)*(kpt(3)+dble(i3))**2&
           & +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
           & +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
           & +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))

    end function dsq
  end subroutine nDimsFromEcut

end module esl_grid_m
