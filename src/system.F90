module system_esl
  use parser_esl

  implicit none
  private

  public ::                &
            system_t,      &
            system_init,   &
            system_end
  
  !Data structure for the system
  type system_t
    !Atoms

    !Cell

  end type system_t

  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine system_init(this, parser)
     type(system_t), intent(inout) :: this
     type(parser_t),  intent(in) :: parser

     !Parse the informations from the input file

   end subroutine system_init
 
   !Release
   !----------------------------------------------------
   subroutine system_end(this)
     type(system_t), intent(inout) :: this

   end subroutine system_end

end module system_esl
