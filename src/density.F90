module density_esl
  use prec, only : dp,ip

  use basis_esl
  use grid_esl

  implicit none
  private

  public :: density_t
  
  !Data structure for the density
  type density_t
    real(kind=dp), allocatable :: density(:)
    real(kind=dp), allocatable :: density_matrix(:,:) 
   contains
     procedure, public :: init
     procedure, public :: calculate
     final  :: cleanup
  end type density_t

  contains

   !Initialize the density
   !----------------------------------------------------
   subroutine init(this, basis, grid)
     class(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis
     type(grid_t),    intent(in)    :: grid

     !Parse the informations from the input file
     select case(basis%basis_type)
       case(PLANEWAVES)
       !Initialization structures for the PW case
       case(ATOMICORBS)
       !Initialization structures for the LO case
     end select 

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(density_t), intent(inout) :: this

     if(allocated(this%density_matrix)) deallocate(this%density_matrix)
     if(allocated(this%density)) deallocate(this%density)

   end subroutine cleanup

   !Calc density
   !----------------------------------------------------
   subroutine calculate(this, basis)
     class(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis

     select case(basis%basis_type)
       case(PLANEWAVES)

       case(ATOMICORBS)

     end select

   end subroutine calculate

end module density_esl
