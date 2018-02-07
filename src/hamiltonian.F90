module hamiltonian_esl
 use parser_esl

 implicit none
 private

 public ::                   &
           hamiltonian_t,    &
           hamiltonian_init, &
           hamiltonian_end,  &
           hamiltonian_apply
          
 !Data structure for the Hamiltonian
 type hamiltonian_t

 end type hamitonian_t

 contains,

   !Initialize the Hamiltonian
   !----------------------------------------------------
   subroutine hamiltonian_init(this, parser)
     type(hamiltonian_t), intent(inout) :: this
     type(parser_t),  intent(in) :: parser

   end subroutine hamiltonian_init


   !Release the Hamiltonian
   !----------------------------------------------------
   subroutine hamiltonian_end(this)
     type(hamiltonian_t), intent(inout) :: this

   end subroutine hamiltonian_end


   !Apply the Hamiltonian matrix
   !----------------------------------------------------
   subroutine hamiltonian_apply(this)
     type(hamiltonian_t), intent(in) :: this

   end subroutine hamiltonian_apply

end module hamiltonian_esl
