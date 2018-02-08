module scf_esl
  use prec, only : dp,ip
  use fdf, only : fdf_get

  use basis_esl
  use density_esl
  use hamiltonian_esl
  use potential_esl
  use system_esl

  implicit none
  private

  public ::                       &
            scf_t,                &
            scf_loop

  !Data structure containing the data for the SCF
  type scf_t 
   real(kind=dp) :: tol_reldens 
   integer(kind=ip) :: max_iter !< Maximum number of iterations
   contains
     private
     procedure, public :: init
     final :: cleanup
  end type scf_t

  contains
 

 !Initialize density and wfn
 !----------------------------------------------------
 subroutine init(this)
   class(scf_t)  :: this

   this%tol_reldens = fdf_get('SCFTolerance',1.0e-6_dp)
   this%max_iter    = fdf_get('SCFMaxIterations', 100)
   !Parse here the data for the SCF
 end subroutine init

 !Cleaning up
 !----------------------------------------------------
 subroutine cleanup(this)
   type(scf_t) :: this


 end subroutine cleanup

 !Perform the self-consistent field calculation
 !----------------------------------------------------
 subroutine scf_loop(this, elsi, hamiltonian, system, states, smear)
  use smear_esl
  use states_esl
    use elsi_wrapper_esl
    type(scf_t),         intent(inout) :: this
    type(elsi_t), intent(inout) :: elsi
   type(hamiltonian_t), intent(inout) :: hamiltonian
   type(system_t),         intent(in) :: system
   type(states_t),         intent(in) :: states
   type(smear_t), intent(inout) :: smear
  
   integer :: iter !< Interation

   do iter = 1, this%max_iter
     !Diagonalization (ELSI/KSsolver)
 
     !Update occupations
     call smear_calc_fermi_and_occ(smear, elsi, states)

     !Calc. density
     call hamiltonian%density%calculate(system%basis)

     !Calc. potentials
     call hamiltonian%potentials%calculate(hamiltonian%density, hamiltonian%energy)

     !Calc. energies
     call hamiltonian%energy%calculate()  

     !Test tolerance and print status

     !Mixing (BLAS/LAPACK)

     !Update Hamiltonian matrix

   end do

 end subroutine scf_loop

end module scf_esl
