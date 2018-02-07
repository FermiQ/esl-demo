module parser_esl

 implicit none
 private

 public ::                   &
           parser_t,         &
           parser_init,      &
           parser_end
 
 !The data structure for the parser         
 type parser_t

 end type parser_t

 contains

   !Initialize the parser
   !----------------------------------------------------
   subroutine parser_init(this)
     type(parser_t), intent(inout) :: parser

   end subroutine parser_init

   !Release the parser
   !----------------------------------------------------
   subroutine parser_end(this)
     type(parser_t), intent(inout) :: parser

   end subroutine parser_end


end module parser_esl
