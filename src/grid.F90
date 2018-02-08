module grid_esl
  use prec, only : dp,ip

  use basis_esl

 implicit none
 private

 public ::  grid_t

 !Data structure for the real space grid
 type grid_t
   real(kind=dp) :: hgrid(3)
   integer :: ndims(3) 
   contains
    private
    procedure, public :: init
    final  :: cleanup
 end type grid_t

 contains

   !Initialize the grid
   !----------------------------------------------------
   subroutine init(this, basis, acell)
     class(grid_t) :: this
     type(basis_t), intent(in) :: basis
     real(kind=dp), intent(in) :: acell
 
     integer :: idim

     !For the moment the spacing in real space is hardcoded
     !For planewave, this must come from the number of G vectors
     this%hgrid(1:3) = 0.25
     do idim = 1,3
       this%ndims(idim) = acell/this%hgrid(idim)
     end do

   end subroutine init


   !Release the grid
   !----------------------------------------------------
   subroutine cleanup(this)
     type(grid_t) :: this

   end subroutine cleanup


end module grid_esl
