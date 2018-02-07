module density_esl
  use prec, only : dp,ip

  use basis_esl

  implicit none
  private

  public ::                 &
            density_t,      &
            density_init,   &
            density_calc,   &
            density_end
  
  !Data structure for the density
  type density_t
    real(kind=dp), allocatable :: density(:)
    real(kind=dp), allocatable :: density_matrix(:,:) 
  end type density_t

  contains

   !Initialize the density
   !----------------------------------------------------
   subroutine density_init(this, basis)
     type(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis

     !Parse the informations from the input file
     select case(basis%basis_type)
       case(PLANEWAVES)
       !Initialization structures for the PW case
       case(ATOMICORBS)
       !Initialization structures for the LO case
     end select 

   end subroutine density_init
 
   !Release
   !----------------------------------------------------
   subroutine density_end(this)
     type(density_t), intent(inout) :: this

     if(allocated(this%density_matrix)) deallocate(this%density_matrix)
     if(allocated(this%density)) deallocate(this%density)

   end subroutine density_end

   !Calc density
   !----------------------------------------------------
   subroutine density_calc(this, basis)
     type(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis

     select case(basis%basis_type)
       case(PLANEWAVES)

       case(ATOMICORBS)

     end select

   end subroutine density_calc 

end module density_esl
