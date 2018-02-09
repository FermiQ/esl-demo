module hamiltonian_esl
  use prec, only : dp,ip

  use density_esl
  use energy_esl
  use potential_esl
  use states_esl
  use system_esl

  implicit none
  private

  public ::                          &
       hamiltonian_t,                &
       hamiltonian_apply,            &
       hamiltonian_apply_local

  !Data structure for the Hamiltonian
  type hamiltonian_t
     type(density_t)   :: density
     type(energy_t)    :: energy
     type(potential_t) :: potentials
   contains
     private
     procedure, public :: init
  end type hamiltonian_t

  interface hamiltonian_apply
    module procedure hamiltonian_dapply, hamiltonian_zapply
  end interface hamiltonian_apply

  interface hamiltonian_apply_local
    module procedure hamiltonian_dapply_local, hamiltonian_zapply_local
  end interface hamiltonian_apply_local

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
   subroutine hamiltonian_dapply(this, psi, hpsi)
     type(hamiltonian_t), intent(in)    :: this
     real(kind=dp),       intent(in)    :: psi(:)
     real(kind=dp),       intent(inout) :: hpsi(:)

     !TODO: Here perform FFT-1
     call hamiltonian_dapply_local(this, psi,hpsi)
     !TODO: Here perform FFT

   end subroutine hamiltonian_dapply

   !Apply the Hamiltonian matrix
   !----------------------------------------------------
   subroutine hamiltonian_zapply(this, psi, hpsi)
     type(hamiltonian_t),   intent(in)    :: this
     complex(kind=dp),      intent(in)    :: psi(:)
     complex(kind=dp),      intent(inout) :: hpsi(:)

     !TODO: Here perform FFT-1
     call hamiltonian_zapply_local(this, psi, hpsi)
     !TODO: Here perform FFT

   end subroutine hamiltonian_zapply

   !Apply the local part of the Hamitonian to a wavefunction
   !----------------------------------------------------
   subroutine hamiltonian_dapply_local(this, psi, hpsi)
     type(hamiltonian_t),  intent(in)    :: this
     real(kind=dp),        intent(in)    :: psi(:)
     real(kind=dp),        intent(inout) :: hpsi(:)

     integer :: ip

     !TODO: Here there is no spin
     forall(ip = 1:this%potentials%np)
       hpsi(ip) = hpsi(ip) + (this%potentials%external(ip) + this%potentials%hartree(ip) + this%potentials%vxc(ip,1))*psi(ip)
     end forall 

   end subroutine hamiltonian_dapply_local


   !Apply the local part of the Hamitonian to a wavefunction
   !----------------------------------------------------
   subroutine hamiltonian_zapply_local(this, psi, hpsi)
     type(hamiltonian_t),   intent(in)    :: this
     complex(kind=dp),      intent(in)    :: psi(:)
     complex(kind=dp),      intent(inout) :: hpsi(:)

     integer :: ip

     !TODO: Here there is no spin
     forall(ip = 1:this%potentials%np)
       hpsi(ip) = hpsi(ip) + (this%potentials%external(ip) + this%potentials%hartree(ip) + this%potentials%vxc(ip,1))*psi(ip)
     end forall

   end subroutine hamiltonian_zapply_local


end module hamiltonian_esl
