module esl_hamiltonian_m
  use prec, only : dp,ip
  use esl_basis_m
  use esl_density_m
  use esl_energy_m
  use esl_force_m
  use esl_geometry_m
  use esl_grid_m
  use esl_ion_interaction_m
  use esl_potential_m
  use esl_states_m

  implicit none
  private

  public ::                          &
      hamiltonian_t,                &
      hamiltonian_apply,            &
      hamiltonian_apply_local

  !Data structure for the Hamiltonian
  type hamiltonian_t
    type(density_t)         :: density
    type(energy_t)          :: energy
    type(force_t)           :: force
    type(ion_interaction_t) :: ion_inter
    type(potential_t)       :: potentials
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
  subroutine init(this, grid, geo, states, periodic)
    class(hamiltonian_t) :: this
    type(grid_t),     intent(in) :: grid
    type(geometry_t), intent(in) :: geo
    type(states_t),   intent(in) :: states
    logical,          intent(in) :: periodic

    call this%density%init(grid)
    call this%energy%init()
    call this%potentials%init(grid, states, periodic)
    call this%force%init(geo%n_atoms)
    call this%ion_inter%init()

    if (periodic) then
      call this%ion_inter%calculate_periodic(geo, this%force%ionion, this%energy%ionion)
    else
      call this%ion_inter%calculate_isolated(geo, this%force%ionion, this%energy%ionion)
    end if

  end subroutine init

  !Apply the Hamiltonian matrix
  !----------------------------------------------------
  subroutine hamiltonian_dapply(this, psi, hpsi)
    type(hamiltonian_t), intent(in)    :: this
    real(dp),       intent(in)    :: psi(:)
    real(dp),       intent(inout) :: hpsi(:)

    !TODO: Here perform FFT-1
    call hamiltonian_dapply_local(this, psi,hpsi)
    !TODO: Here perform FFT

  end subroutine hamiltonian_dapply

  !Apply the Hamiltonian matrix
  !----------------------------------------------------
  subroutine hamiltonian_zapply(this, psi, hpsi)
    type(hamiltonian_t),   intent(in)    :: this
    complex(dp),      intent(in)    :: psi(:)
    complex(dp),      intent(inout) :: hpsi(:)

    !TODO: Here perform FFT-1
    call hamiltonian_zapply_local(this, psi, hpsi)
    !TODO: Here perform FFT

  end subroutine hamiltonian_zapply

  !Apply the local part of the Hamitonian to a wavefunction
  !----------------------------------------------------
  subroutine hamiltonian_dapply_local(this, psi, hpsi)
    type(hamiltonian_t),  intent(in)    :: this
    real(dp),        intent(in)    :: psi(:)
    real(dp),        intent(inout) :: hpsi(:)

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
    complex(dp),      intent(in)    :: psi(:)
    complex(dp),      intent(inout) :: hpsi(:)

    integer :: ip

    !TODO: Here there is no spin
    forall(ip = 1:this%potentials%np)
      hpsi(ip) = hpsi(ip) + (this%potentials%external(ip) + this%potentials%hartree(ip) + this%potentials%vxc(ip,1))*psi(ip)
    end forall

  end subroutine hamiltonian_zapply_local

end module esl_hamiltonian_m
