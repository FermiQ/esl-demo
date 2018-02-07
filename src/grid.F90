module grid_esl
  use prec, only : dp,ip

 implicit none
 private

 public ::                 &
           grid_t,         &
           grid_init,      &
           grid_end

 !Data structure for the real space grid
 type grid_t

 end type grid_t

 contains

   !Initialize the grid
   !----------------------------------------------------
   subroutine grid_init(this)
     type(grid_t), intent(inout) :: this

   end subroutine grid_init


   !Release the grid
   !----------------------------------------------------
   subroutine grid_end(this)
     type(grid_t):: this

   end subroutine grid_end

end module grid_esl
