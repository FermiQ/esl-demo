module density_esl

  implicit none
  private

  public ::                 &
            density_t,      &
            density_init,   &
            density_end
  
  !Data structure for the density
  type density_t

  end type density_t

  contains

   !Initialize the density
   !----------------------------------------------------
   subroutine density_init(this)
     type(density_t), intent(inout) :: this

     !Parse the informations from the input file

   end subroutine density_init
 
   !Release
   !----------------------------------------------------
   subroutine density_end(this)
     type(density_t), intent(inout) :: this

   end subroutine density_end

end module density_esl
