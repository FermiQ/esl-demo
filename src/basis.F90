module basis_esl
  use parser_esl

  implicit none
  private

  public ::                &
            basis_t,      &
            basis_init,   &
            basis_end
  
  !Data structure for the basis
  type basis_t
    integer :: basis_type
  end type basis_t

  integer, public, parameter :: &
    PLANEWAVES   = 1,           &
    ATOMICORBS   = 2


  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine basis_init(this, parser)
     type(basis_t), intent(inout) :: this
     type(parser_t),  intent(in) :: parser

     !Parse the informations from the input file

   end subroutine basis_init
 
   !Release
   !----------------------------------------------------
   subroutine basis_end(this)
     type(basis_t), intent(inout) :: this

   end subroutine basis_end

end module basis_esl
