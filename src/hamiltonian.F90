module hamiltonian_esl

 implicit none
 private

 public ::                   &
           hamiltonian_t,    &
           hamiltonian_apply
          
 !Data structure for the Hamiltonian
 type hamiltonian_t
    integer :: dummy ! remove once proper data members are added
    contains
      private
      procedure, public :: init
      final :: cleanup
 end type 

 contains

   !Initialize the Hamiltonian
   !----------------------------------------------------
   subroutine init(this)
     class(hamiltonian_t) :: this

   end subroutine init


   !Release the Hamiltonian
   !----------------------------------------------------
   subroutine cleanup(this)
     type(hamiltonian_t):: this

   end subroutine cleanup


   !Apply the Hamiltonian matrix
   !----------------------------------------------------
   subroutine hamiltonian_apply(this)
     type(hamiltonian_t), intent(in) :: this

   end subroutine hamiltonian_apply

end module hamiltonian_esl
