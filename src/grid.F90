module grid_esl
  use prec, only : dp,ip

  use basis_esl

 implicit none
 private

 public ::  grid_t

 !Data structure for the real space grid
 type grid_t
   real(kind=dp) :: hgrid(3) !< Real space spacing
   integer :: ndims(3)  !< Number of points in each directions
   integer :: np !< Total number of points in the real space grid
   contains
    private
    procedure, public :: init
    procedure, public :: summary
    final  :: cleanup
 end type grid_t

 contains

   !Initialize the grid
   !----------------------------------------------------
   subroutine init(this, basis, cell)
     use module_fft_sg
     class(grid_t) :: this
     type(basis_t), intent(in) :: basis
     real(kind=dp), dimension(3,3), intent(in) :: cell
 
     integer :: idim, n, twice

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


end module grid_esl
