module hamiltonian_esl
  use density_esl
  use basis_esl

 implicit none
 private

 public ::                   &
           hamiltonian_t,    &
           hamiltonian_apply
          
 !Data structure for the Hamiltonian
 type hamiltonian_t
   type(density_t) :: density
    contains
      private
      procedure, public :: init
      final :: cleanup
 end type 

 contains

   !Initialize the Hamiltonian
   !----------------------------------------------------
   subroutine init(this, basis)
     class(hamiltonian_t) :: this
     type(basis_t), intent(in) :: basis

     call density_init(this%density, basis)

   end subroutine init


   !Release the Hamiltonian
   !----------------------------------------------------
   subroutine cleanup(this)
     type(hamiltonian_t):: this

     call density_end(this%density)

   end subroutine cleanup


   !Apply the Hamiltonian matrix
   !----------------------------------------------------
   subroutine hamiltonian_apply(this)
     type(hamiltonian_t), intent(in) :: this

   end subroutine hamiltonian_apply

end module hamiltonian_esl
