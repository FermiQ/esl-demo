module grid_esl
  use prec, only : dp,ip

  use basis_esl
  use numeric_esl, only : grylmr

 implicit none
 private

 public ::           &
           grid_t,   &
           integrate

 !Data structure for the real space grid
 type grid_t
   real(kind=dp) :: hgrid(3) !< Real space spacing
   integer :: ndims(3)  !< Number of points in each directions
   integer :: np !< Total number of points in the real space grid
   real(kind=dp), allocatable :: r(:,:) !<Grid point coordinates 
   real(kind=dp) :: volelem !<Volume element
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
   subroutine init(this, basis, cell)
     use module_fft_sg
     class(grid_t) :: this
     type(basis_t), intent(in) :: basis
     real(kind=dp), dimension(3,3), intent(in) :: cell
 
     integer :: idim, ix, iy, iz, ip
     integer :: n, twice

     !For the moment the spacing in real space is hardcoded
     !For planewave, this must come from the number of G vectors
     this%hgrid(1:3) = .25
     do idim = 1,3
        n = ceiling(cell(idim, idim) / this%hgrid(idim)) - 1
        do
           call fourier_dim(n, this%ndims(idim))
           call fourier_dim(2 * n, twice)
           this%hgrid(idim) = cell(idim, idim) / real(this%ndims(idim) + 1, dp)
           if (2 * n == twice) exit ! Ensure that Poisson Solver double grid will work.
        end do
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
     real(kind=dp), intent(out) :: r_at(3)
     real(kind=dp), intent(out) :: ao(:)
     real(kind=dp), intent(out) :: grad_ao(:,:)

     integer :: ip
     real(kind=dp) :: x, y, z, r

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
     real(kind=dp),   intent(in) :: ff(:)
     real(kind=dp),  intent(out) :: int_ff

     integer :: ip

     int_ff = 0.d0
     forall(ip=1:grid%np)
       int_ff = int_ff + ff(ip)
     end forall
     int_ff = int_ff*grid%volelem

   end subroutine dintegrate

   !Integrate a function over the real-space grid
   !----------------------------------------------------
   subroutine zintegrate(grid, ff, int_ff)
     type(grid_t),       intent(in) :: grid
     complex(kind=dp),   intent(in) :: ff(:)
     complex(kind=dp),  intent(out) :: int_ff

     integer :: ip

     int_ff = cmplx(0.d0,0.d0)
     forall(ip=1:grid%np)
       int_ff = int_ff + ff(ip)
     end forall
     int_ff = int_ff*grid%volelem

   end subroutine zintegrate


end module grid_esl
