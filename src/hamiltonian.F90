module hamiltonian_esl
  use density_esl
  use energy_esl
  use potential_esl
  use system_esl

 implicit none
 private

 public ::                   &
           hamiltonian_t,    &
           hamiltonian_apply
          
 !Data structure for the Hamiltonian
 type hamiltonian_t
   type(density_t)   :: density
   type(energy_t)    :: energy
   type(potential_t) :: potentials
    contains
      private
      procedure, public :: init
      final :: cleanup
 end type 

 contains

   !Initialize the Hamiltonian
   !----------------------------------------------------
   subroutine init(this, sys)
     class(hamiltonian_t) :: this
     type(system_t), intent(in) :: sys

     call density_init(this%density, sys%basis, sys%grid)
     call this%energy%init()
     call potential_init(this%potentials, sys%basis, sys%grid)

   end subroutine init


   !Release the Hamiltonian
   !----------------------------------------------------
   subroutine cleanup(this)
     type(hamiltonian_t):: this

     call density_end(this%density)
     call potential_end(this%potentials)

   end subroutine cleanup


   !Apply the Hamiltonian matrix
   !----------------------------------------------------
   subroutine hamiltonian_apply(this)
     type(hamiltonian_t), intent(in) :: this

   end subroutine hamiltonian_apply

end module hamiltonian_esl
