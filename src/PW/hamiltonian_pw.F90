module esl_hamiltonian_pw_m
  use prec, only : dp,ip

  use esl_potential_m
  use esl_states_m

  implicit none
  private

  public ::                          &
      hamiltonian_pw_t,                &
      hamiltonian_pw_apply,            &
      hamiltonian_pw_apply_local

  !Data structure for the Hamiltonian
  type hamiltonian_pw_t
  contains
    private
    procedure, public :: init
    procedure, public :: eigensolver
    final :: cleanup
  end type hamiltonian_pw_t

  interface hamiltonian_pw_apply
    module procedure hamiltonian_pw_dapply, hamiltonian_pw_zapply
  end interface hamiltonian_pw_apply

  interface hamiltonian_pw_apply_local
    module procedure hamiltonian_pw_dapply_local, hamiltonian_pw_zapply_local
  end interface hamiltonian_pw_apply_local

contains

  !Initialize the Hamiltonian
  !----------------------------------------------------
  subroutine init(this)
    class(hamiltonian_pw_t) :: this

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(hamiltonian_pw_t) :: this

  end subroutine cleanup

  !Eigensolver
  !----------------------------------------------------
  subroutine eigensolver(this, states)
    class(hamiltonian_pw_t) :: this
    type(states_t), intent(inout) :: states

  end subroutine eigensolver


  !Apply the Hamiltonian matrix
  !----------------------------------------------------
  subroutine hamiltonian_pw_dapply(pot, psi, hpsi)
    type(potential_t),  intent(in)    :: pot
    real(dp),       intent(in)        :: psi(:)
    real(dp),       intent(inout)     :: hpsi(:)

    !TODO: Here perform FFT-1
    call hamiltonian_pw_dapply_local(pot, psi,hpsi)
    !TODO: Here perform FFT

  end subroutine hamiltonian_pw_dapply

  !Apply the Hamiltonian matrix
  !----------------------------------------------------
  subroutine hamiltonian_pw_zapply(pot, psi, hpsi)
    type(potential_t),  intent(in)    :: pot
    complex(dp),      intent(in)      :: psi(:)
    complex(dp),      intent(inout)   :: hpsi(:)

    !TODO: Here perform FFT-1
    call hamiltonian_pw_zapply_local(pot, psi, hpsi)
    !TODO: Here perform FFT

  end subroutine hamiltonian_pw_zapply

  !Apply the local part of the Hamitonian to a wavefunction
  !----------------------------------------------------
  subroutine hamiltonian_pw_dapply_local(pot, psi, hpsi)
    type(potential_t),  intent(in)    :: pot
    real(dp),           intent(in)    :: psi(:)
    real(dp),           intent(inout) :: hpsi(:)

    integer :: ip

    !TODO: Here there is no spin
    forall(ip = 1:pot%np)
      hpsi(ip) = hpsi(ip) + (pot%external(ip) + pot%hartree(ip) + pot%vxc(ip))*psi(ip)
    end forall

  end subroutine hamiltonian_pw_dapply_local


  !Apply the local part of the Hamitonian to a wavefunction
  !----------------------------------------------------
  subroutine hamiltonian_pw_zapply_local(pot, psi, hpsi)
    type(potential_t),   intent(in)    :: pot
    complex(dp),            intent(in)    :: psi(:)
    complex(dp),            intent(inout) :: hpsi(:)

    integer :: ip

    !TODO: Here there is no spin
    forall(ip = 1:pot%np)
      hpsi(ip) = hpsi(ip) + (pot%external(ip) + pot%hartree(ip) + pot%vxc(ip))*psi(ip)
    end forall

  end subroutine hamiltonian_pw_zapply_local

end module esl_hamiltonian_pw_m
