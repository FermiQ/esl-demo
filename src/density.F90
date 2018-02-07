module density_esl
  use prec, only : dp,ip

  use basis_esl

  implicit none
  private

  public ::                 &
            density_t,      &
            density_init,   &
            densitu_calc,   &
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
     select(basis%basis_type)
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

     if(allocated(density_matrix)) deallocate(density_matrix)
     if(allocated(density)) deallocate(density)

   end subroutine density_end

   !Calc density
   !----------------------------------------------------
   subroutine density_calc(this, basis)
     type(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis

     select(basis%basis_type)
       case(PLANEWAVES)

       case(ATOMICORBS)

     end select

   end subroutine density_calc 

end module density_esl
