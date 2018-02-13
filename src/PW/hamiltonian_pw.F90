module esl_hamiltonian_pw_m
  use prec, only : dp,ip
  use elsi
  use elsi_rci
  use elsi_rci_constants
  use elsi_rci_omm
  use elsi_rci_precision, only:r8, i4

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
    type(elsi_handle)    :: e_h
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
  subroutine init(this, states)
    class(hamiltonian_pw_t) :: this
    type(states_t), intent(in) :: states

    integer(kind=i4) :: solver, matrix_size

    call elsi_init(this%e_h, solver, 1, 0, matrix_size, states%nel, states%nstates)

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(hamiltonian_pw_t) :: this

    ! Finalize ELSI
    call elsi_finalize(this%e_h)

  end subroutine cleanup

  !Eigensolver
  !----------------------------------------------------
  subroutine eigensolver(this, states)
    class(hamiltonian_pw_t) :: this
    type(states_t), intent(inout) :: states
   
    real(kind=r8), allocatable :: result_in(:)
    type(rci_instr)      :: iS
    integer(kind=i4) :: task, m, n, ijob
    real(kind=r8) :: e_min = 0.0_r8
    real(kind=r8) :: cg_tol = 0.0001_r8
    integer(kind=i4) :: max_iter = 100
    logical :: long_out = .false.

    allocate (result_in(1))
    ijob = -1
    do
        call rci_omm(ijob, iS, task, result_in, m, n, &
                     e_min, cg_tol, max_iter, long_out)

        select case (task)
        case (ELSI_RCI_NULL)
        case (ELSI_RCI_CONVERGE)
            exit
        case (ELSI_RCI_H_MULTI)
   !         call mm_multiply(H, iS%TrH, Work(iS%Aidx), &
   !                          'n', Work(iS%Bidx), 1.0_r8, 0.0_r8, &
   !                          m_operation)
        case (ELSI_RCI_S_MULTI)
   !         call mm_multiply(S, iS%TrS, Work(iS%Aidx), &
   !                          'n', Work(iS%Bidx), 1.0_r8, 0.0_r8, &
   !                          m_operation)
        case (ELSI_RCI_P_MULTI)
            ! No preconditioner
   !         call m_add(Work(iS%Aidx), 'n', Work(iS%Bidx), 1.0_r8, &
   !                    0.0_r8, m_operation)
        case (ELSI_RCI_GEMM)
   !         call mm_multiply(Work(iS%Aidx), iS%TrA, Work(iS%Bidx), &
   !                          iS%TrB, Work(iS%Cidx), iS%alpha, iS%beta, &
   !                          m_operation)
        case (ELSI_RCI_AXPY)
   !         call m_add(Work(iS%Aidx), 'n', Work(iS%Bidx), iS%alpha, &
   !                    1.0_r8, m_operation)
        case (ELSI_RCI_COPY)
   !         call m_add(Work(iS%Aidx), iS%TrA, Work(iS%Bidx), 1.0_r8, &
   !                    0.0_r8, m_operation)
        case (ELSI_RCI_TRACE)
   !         call m_trace(Work(iS%Aidx), result_in(1), m_operation)
        case (ELSI_RCI_DOT)
   !         call mm_trace(Work(iS%Aidx), Work(iS%Bidx), result_in(1), &
   !                       m_operation)
        case (ELSI_RCI_SCALE)
   !         call m_scale(Work(iS%Aidx), iS%alpha, m_operation)
        case default
        end select

    end do

    deallocate (result_in)

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

    !Note that hartree contains the external potential (from PSolver)
    forall(ip = 1:pot%np)
      hpsi(ip) = hpsi(ip) + (pot%hartree(ip) + pot%vxc(ip))*psi(ip)
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
 
    !Note that hartree contains the external potential (from PSolver)
    forall(ip = 1:pot%np)
      hpsi(ip) = hpsi(ip) + (pot%hartree(ip) + pot%vxc(ip))*psi(ip)
    end forall

  end subroutine hamiltonian_pw_zapply_local

end module esl_hamiltonian_pw_m
