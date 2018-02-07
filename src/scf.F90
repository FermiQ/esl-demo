module scf_esl
  use prec, only : dp,ip
  use fdf, only : fdf_integer, fdf_double

  use basis_esl
  use density_esl
  use hamiltonian_esl
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

     this%tol_reldens=fdf_double('SCFTolerance',1.0e-6_dp)
     this%max_iter=fdf_integer('SCFMaxIterations', 100)
   !Parse here the data for the SCF
 end subroutine init

 !Cleaning up
 !----------------------------------------------------
 subroutine cleanup(this)
   type(scf_t) :: this


 end subroutine cleanup

 !Perform the self-consistent field calculation
 !----------------------------------------------------
 subroutine scf_loop(this, hamiltonian, system)
   type(scf_t),         intent(inout) :: this
   type(hamiltonian_t), intent(inout) :: hamiltonian
   type(system_t),         intent(in) :: system
  
   integer :: iter !< Interation

   do iter = 1, this%max_iter
     !Diagonalization (ELSI/KSsolver)
 
     !Update occupations

     !Calc. density
     call density_calc(hamiltonian%density, system%basis)    

     !Test tolerance and print status

     !Mixing (BLAS/LAPACK)

     !Update Hamiltonian matrix

   end do

 end subroutine scf_loop

end module scf_esl
