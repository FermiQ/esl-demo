module grid_esl
  use prec, only : dp,ip

  use basis_esl

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
   real(kind=dp) :: volelem !<Volume element
   contains
    private
    procedure, public :: init
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
 
     integer :: idim, n, twice

     select case (basis%basis_type)
     case (PLANEWAVES)
        call nDimsFromEcut(this%ndims, basis%pw_basis%ecut, cell, [0._dp, 0._dp, 0._dp])
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

     !We have a cubic cell
     this%volelem = this%hgrid(1)*this%hgrid(2)*this%hgrid(3)

   end subroutine init


   !Release the grid
   !----------------------------------------------------
   subroutine cleanup(this)
     type(grid_t) :: this

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

   subroutine nDimsFromEcut(ndims, ecut, cell, kpt)
     use numerics, only: pi
     use yaml_output
     implicit none
     integer, dimension(3), intent(out) :: ndims
     real(dp), intent(in) :: ecut
     real(dp), dimension(3,3), intent(in) :: cell
     real(dp), dimension(3) :: kpt

     real(dp) :: threshold
     integer :: dir, n, i
     real(dp), dimension(3,3) :: gmet, gcell
     real(dp), parameter :: boxcutmin = 2.0_dp

     call matr3inv(cell, gcell)
     do i = 1, 3
        gmet(i, :) = gcell(1, i) * gcell(1, :) + &
             &   gcell(2, i) * gcell(2, :) + &
             &   gcell(3, i) * gcell(3, :)
     end do

     threshold = 0.5_dp * boxcutmin**2 * ecut / pi**2
     !@todo: Don't take into account symmetries or k points.
     ndims = 16
     do
        if (smallest(ndims, gmet, dir) >= threshold) exit
        call fourier_dim(ndims(dir) + 1, ndims(dir))
     end do

   contains
     subroutine matr3inv(matin, matout)
       implicit none
       real(dp), dimension(3,3), intent(in) :: matin
       real(dp), dimension(3,3), intent(out) :: matout

       real(dp) :: dd,t1,t2,t3

       t1 = matin(2,2) * matin(3,3) - matin(3,2) * matin(2,3)
       t2 = matin(3,2) * matin(1,3) - matin(1,2) * matin(3,3)
       t3 = matin(1,2) * matin(2,3) - matin(2,2) * matin(1,3)
       dd  = 1.d0/ (matin(1,1) * t1 + matin(2,1) * t2 + matin(3,1) * t3)
       matout(1,1) = t1 * dd
       matout(2,1) = t2 * dd
       matout(3,1) = t3 * dd
       matout(1,2) = (matin(3,1)*matin(2,3)-matin(2,1)*matin(3,3)) * dd
       matout(2,2) = (matin(1,1)*matin(3,3)-matin(3,1)*matin(1,3)) * dd
       matout(3,2) = (matin(2,1)*matin(1,3)-matin(1,1)*matin(2,3)) * dd
       matout(1,3) = (matin(2,1)*matin(3,2)-matin(3,1)*matin(2,2)) * dd
       matout(2,3) = (matin(3,1)*matin(1,2)-matin(1,1)*matin(3,2)) * dd
       matout(3,3) = (matin(1,1)*matin(2,2)-matin(2,1)*matin(1,2)) * dd
     end subroutine matr3inv

     function smallest(ndims, gmet, dir)
       implicit none
       integer, dimension(3), intent(in) :: ndims
       real(dp), dimension(3,3), intent(in) :: gmet
       integer, intent(out) :: dir
       real(dp) :: smallest

       integer :: i1, i2, i3, alpha, beta, gamma, s1, s2, s3, idir
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
       implicit none
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
end module grid_esl
