module hamiltonian_esl
  use density_esl
  use energy_esl
  use potential_esl
  use states_esl
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
  end type hamiltonian_t

contains

  !Initialize the Hamiltonian
  !----------------------------------------------------
   subroutine init(this, sys, states)
     class(hamiltonian_t) :: this
     type(system_t), intent(in) :: sys
     type(states_t), intent(in) :: states

     call this%density%init(sys%basis, sys%grid)
     call this%energy%init()
     call this%potentials%init(sys%basis, sys%grid, states)

   end subroutine init

   !Apply the Hamiltonian matrix
   !----------------------------------------------------
   subroutine hamiltonian_apply(this)
     type(hamiltonian_t), intent(in) :: this

   end subroutine hamiltonian_apply

end module hamiltonian_esl
